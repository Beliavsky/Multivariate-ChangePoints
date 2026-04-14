/*
 * xreturns_corr_cp.cpp
 *
 * C++ equivalent of xreturns_corr_cp.f90.
 * Reads asset prices, computes percent returns, and for each asset pair
 * finds changepoints in correlation using O(n^2 * max_m) dynamic programming.
 *
 * Cost matrix layout: cost[end*n + start] (end-major).
 * This gives sequential memory access in both construction (outer j=end,
 * inner i=start) and the DP inner loop (fixed end=i, varying start=k+1).
 *
 * Compile: g++ -O2 -std=c++17 -o xreturns_corr_cp xreturns_corr_cp.cpp
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ── parameters ────────────────────────────────────────────────────────────────
static const std::string              PRICES_FILE = "spy_efa_eem_tlt.csv";
static const double                   SCALE_RET   = 100.0;
static const int                      MAX_CP      = 20;
static const int                      MIN_SEG_LEN = 50;
static const bool                     PRINT_PARAMS = false;
// Empty → use all columns; otherwise restrict to listed tickers.
static const std::vector<std::string> SYM_ALLOWED = {};  // e.g. {"SPY","TLT"}
// ─────────────────────────────────────────────────────────────────────────────

static constexpr double BIG = 1.0e20;

// ── CSV reader ────────────────────────────────────────────────────────────────

static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ','))
        out.push_back(tok);
    return out;
}

struct CsvData {
    std::vector<std::string>         dates;    // one entry per price row
    std::vector<std::string>         columns;  // asset names
    std::vector<std::vector<double>> prices;   // prices[col][row]
};

static CsvData read_csv(const std::string& filename,
                        const std::vector<std::string>& sym_allowed) {
    std::ifstream f(filename);
    if (!f) { std::cerr << "Cannot open " << filename << "\n"; std::exit(1); }

    CsvData data;
    std::string line;

    std::getline(f, line);
    auto hdr = split_csv(line);

    std::vector<int> col_idx;
    for (int c = 1; c < (int)hdr.size(); ++c) {
        bool keep = sym_allowed.empty() ||
            std::find(sym_allowed.begin(), sym_allowed.end(), hdr[c]) != sym_allowed.end();
        if (keep) {
            col_idx.push_back(c);
            data.columns.push_back(hdr[c]);
        }
    }
    data.prices.resize(col_idx.size());

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto fields = split_csv(line);
        data.dates.push_back(fields[0]);
        for (int c = 0; c < (int)col_idx.size(); ++c) {
            int idx = col_idx[c];
            data.prices[c].push_back(
                idx < (int)fields.size() ? std::stod(fields[idx]) : 0.0);
        }
    }
    return data;
}

// ── returns ───────────────────────────────────────────────────────────────────

static void compute_returns(const CsvData& csv, double scale,
                             std::vector<std::vector<double>>& ret,
                             std::vector<std::string>& ret_dates) {
    int ncols = (int)csv.prices.size();
    int nrows = (int)csv.prices[0].size();
    int n     = nrows - 1;

    ret.assign(ncols, std::vector<double>(n));
    for (int c = 0; c < ncols; ++c)
        for (int t = 0; t < n; ++t)
            ret[c][t] = scale * (csv.prices[c][t+1] / csv.prices[c][t] - 1.0);

    ret_dates.assign(csv.dates.begin() + 1, csv.dates.end());
}

// ── cost matrix ───────────────────────────────────────────────────────────────

// cost[j*n + i] = 0.5*(j-i+1)*log(1-r²) for segment i..j (0-indexed, inclusive).
// Outer loop over j=end, inner over i=start → sequential writes within each row j.
// BIG for segments shorter than min_seg_len or where j < i.
static void fill_cost_matrix(std::vector<double>& cost,
                              const std::vector<double>& x,
                              const std::vector<double>& y,
                              int n, int min_seg_len) {
    std::fill(cost.begin(), cost.end(), BIG);

    std::vector<double> sx(n+1,0), sy(n+1,0), sxx(n+1,0), syy(n+1,0), sxy(n+1,0);
    for (int t = 0; t < n; ++t) {
        sx[t+1]  = sx[t]  + x[t];
        sy[t+1]  = sy[t]  + y[t];
        sxx[t+1] = sxx[t] + x[t]*x[t];
        syy[t+1] = syy[t] + y[t]*y[t];
        sxy[t+1] = sxy[t] + x[t]*y[t];
    }

    for (int j = min_seg_len - 1; j < n; ++j) {
        double* row  = &cost[(size_t)j * n];   // sequential writes: row[0], row[1], ...
        double sx_j  = sx[j+1],  sy_j  = sy[j+1];
        double sxx_j = sxx[j+1], syy_j = syy[j+1], sxy_j = sxy[j+1];

        for (int i = 0; i <= j - min_seg_len + 1; ++i) {
            double m     = j - i + 1;
            double sx_s  = sx_j  - sx[i];
            double sy_s  = sy_j  - sy[i];
            double sxx_s = sxx_j - sxx[i];
            double syy_s = syy_j - syy[i];
            double sxy_s = sxy_j - sxy[i];

            double ma = sx_s / m,  mb = sy_s / m;
            double va = sxx_s / m - ma*ma;
            double vb = syy_s / m - mb*mb;
            double r  = (sxy_s / m - ma*mb) / std::sqrt(std::max(va*vb, 1.0e-20));
            r = std::max(-1.0+1e-10, std::min(1.0-1e-10, r));

            row[i] = 0.5 * m * std::log(std::max(1.0 - r*r, 1.0e-10));
        }
    }
}

// ── dynamic programming ───────────────────────────────────────────────────────

// dp[m*n + i]  = min total cost for obs 0..i split into m+1 segments (m 0-indexed).
// par[m*n + i] = 0-indexed last obs k of the m-th predecessor segment.
//
// For fixed i and the inner k-loop: dp_prev[k] and cost[i*n + k+1] are both
// accessed sequentially → good cache behaviour.
static void solve_dp(const std::vector<double>& cost, int n, int max_m,
                     std::vector<double>& dp, std::vector<int>& par) {
    dp.assign((size_t)max_m * n, BIG);
    par.assign((size_t)max_m * n, 0);

    // m=0 (1 segment): segment 0..i → cost[i*n + 0]
    for (int i = 0; i < n; ++i)
        dp[i] = cost[(size_t)i * n];

    for (int m = 1; m < max_m; ++m) {
        double*       dp_cur  = &dp[(size_t)m * n];
        int*          par_cur = &par[(size_t)m * n];
        const double* dp_prev = &dp[(size_t)(m-1) * n];

        for (int i = m; i < n; ++i) {
            // k runs from m-1 to i-1; both arrays accessed sequentially.
            const double* cost_row = &cost[(size_t)i * n + m];  // cost[i][k+1], k=m-1..i-1
            const double* dp_k     = dp_prev + (m-1);           // dp[m-1][k],   k=m-1..i-1
            int nk = i - m + 1;

            double best_val = BIG;
            int    best_ki  = 0;
            for (int ki = 0; ki < nk; ++ki) {
                double val = dp_k[ki] + cost_row[ki];
                if (val < best_val) { best_val = val; best_ki = ki; }
            }
            dp_cur[i]  = best_val;
            par_cur[i] = (m-1) + best_ki;   // k = (m-1) + ki
        }
    }
}

// ── changepoint reconstruction ────────────────────────────────────────────────

static std::vector<int> reconstruct_cps(const std::vector<int>& par,
                                         int n, int m_segs) {
    std::vector<int> cps;
    int cp = n - 1;
    for (int k = m_segs - 1; k >= 1; --k) {
        cp = par[(size_t)k * n + cp];
        cps.push_back(cp);
    }
    std::sort(cps.begin(), cps.end());
    return cps;
}

// ── segment correlation ───────────────────────────────────────────────────────

static double seg_corr(const std::vector<double>& x, const std::vector<double>& y,
                        int start, int end) {
    int m = end - start + 1;
    double mx = 0, my = 0;
    for (int i = start; i <= end; ++i) { mx += x[i]; my += y[i]; }
    mx /= m; my /= m;
    double num = 0, dx2 = 0, dy2 = 0;
    for (int i = start; i <= end; ++i) {
        double dx = x[i]-mx, dy = y[i]-my;
        num += dx*dy; dx2 += dx*dx; dy2 += dy*dy;
    }
    double den = std::sqrt(dx2 * dy2);
    return den > 1e-15 ? num / den : 0.0;
}

// ── output ────────────────────────────────────────────────────────────────────

static void print_pair_results(const std::vector<double>& dp,
                                const std::vector<int>&    par,
                                const std::vector<double>& x,
                                const std::vector<double>& y,
                                int n, int max_m,
                                const std::vector<std::string>& dates,
                                bool print_params) {
    std::cout << std::setw(2)  << "M"
              << std::setw(12) << "LL"
              << std::setw(12) << "AIC"
              << std::setw(12) << "BIC"
              << "    ChangePoints\n";

    double log_n   = std::log(n);
    double min_aic = BIG, min_bic = BIG;
    int    best_aic = 1, best_bic = 1;

    for (int m_segs = 1; m_segs <= max_m; ++m_segs) {
        int    m          = m_segs - 1;
        double total_cost = dp[(size_t)m * n + (n-1)];
        if (total_cost >= 1e19) continue;

        double ll       = -total_cost;
        int    k_params = 2*m_segs - 1;
        double aic      = -2.0*ll + 2.0*k_params;
        double bic      = -2.0*ll + (double)k_params * log_n;

        if (aic < min_aic) { min_aic = aic; best_aic = m_segs; }
        if (bic < min_bic) { min_bic = bic; best_bic = m_segs; }

        auto cps = reconstruct_cps(par, n, m_segs);

        std::cout << std::setw(2) << m_segs
                  << std::fixed << std::setprecision(2)
                  << std::setw(12) << ll
                  << std::setw(12) << aic
                  << std::setw(12) << bic;
        for (int cp : cps)
            std::cout << "  " << (cp+1);  // 1-indexed, matches Fortran output
        std::cout << "\n";

        if (print_params) {
            std::vector<int> bounds = {-1};
            bounds.insert(bounds.end(), cps.begin(), cps.end());
            bounds.push_back(n-1);
            std::cout << "  " << std::setw(12) << "Start"
                      << std::setw(12) << "End"
                      << std::setw(6)  << "#obs"
                      << std::setw(9)  << "Corr\n";
            for (int k = 0; k < m_segs; ++k) {
                int s = bounds[k]+1, e = bounds[k+1];
                std::cout << "  " << std::setw(12) << dates[s]
                          << std::setw(12) << dates[e]
                          << std::setw(6)  << (e-s+1)
                          << std::fixed << std::setprecision(4)
                          << std::setw(9) << seg_corr(x,y,s,e) << "\n";
            }
        }
    }

    std::cout << "Changepoints chosen by AIC, BIC: "
              << (best_aic-1) << "  " << (best_bic-1) << "\n";
}

// ── main ──────────────────────────────────────────────────────────────────────

int main() {
    auto t0 = std::chrono::steady_clock::now();

    auto csv   = read_csv(PRICES_FILE, SYM_ALLOWED);
    int n_days = (int)csv.dates.size();
    int n_cols = (int)csv.columns.size();
    int n      = n_days - 1;

    std::cout << "read " << n_days << " days and " << n_cols
              << " columns from " << PRICES_FILE << "\n"
              << std::fixed << std::setprecision(4)
              << "return scaling: " << SCALE_RET << "\n";

    std::vector<std::vector<double>> ret;
    std::vector<std::string>         ret_dates;
    compute_returns(csv, SCALE_RET, ret, ret_dates);

    int max_m = MAX_CP + 1;

    // Pre-allocate; cost is n×n ~250 MB for n=5587, reused across pairs.
    std::vector<double> cost((size_t)n * n);
    std::vector<double> dp;
    std::vector<int>    par;

    for (int i = 0; i < n_cols; ++i) {
        for (int j = i+1; j < n_cols; ++j) {
            std::cout << "\nchanges in "
                      << csv.columns[i] << "-" << csv.columns[j]
                      << " correlation\n";

            auto t_pair = std::chrono::steady_clock::now();
            fill_cost_matrix(cost, ret[i], ret[j], n, MIN_SEG_LEN);
            solve_dp(cost, n, max_m, dp, par);
            print_pair_results(dp, par, ret[i], ret[j],
                               n, max_m, ret_dates, PRINT_PARAMS);
            auto t1 = std::chrono::steady_clock::now();
            std::cout << "  ("
                      << std::fixed << std::setprecision(1)
                      << std::chrono::duration<double>(t1 - t_pair).count()
                      << "s)\n";
        }
    }

    auto t_end = std::chrono::steady_clock::now();
    std::cout << "\nwall time elapsed: "
              << std::fixed << std::setprecision(2)
              << std::chrono::duration<double>(t_end - t0).count() << "s\n";
    return 0;
}
