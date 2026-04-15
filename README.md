# Multivariate-ChangePoints

Fortran programs (with Python and C++ equivalents for the main program) that find structural breaks in the mean, variance, correlation, and covariance of financial return series using exact dynamic programming.

Two sample data files are included:

**`spy_efa_eem_tlt.csv`** — daily closing prices for four ETFs from 1993 to 2015:

| Ticker | Description |
|--------|-------------|
| SPY | S&P 500 |
| EFA | Developed international equities |
| EEM | Emerging market equities |
| TLT | 20+ year US Treasury bonds |

**`asset_class_etf_prices.csv`** — daily closing prices for nine ETFs from 2007 to 2026:

| Ticker | Description |
|--------|-------------|
| SPY | S&P 500 |
| EFA | Developed international equities |
| EEM | Emerging market equities |
| EMB | Emerging market bonds |
| HYG | High-yield corporate bonds |
| LQD | Investment-grade corporate bonds |
| TLT | 20+ year US Treasury bonds |
| GLD | Gold |
| USO | Oil |

---

## Statistical models

All changepoint programs find the exact global optimum for each fixed number of segments m using the recurrence

    dp(i, m) = min over k < i of [ dp(k, m-1) + cost(k+1, i) ]

in O(n² × max_m) time, then select the number of segments by AIC and BIC.  The cost functions and parameter counts differ by model:

| Program | Series z_t | Cost(i..j) | Params/segment |
|---------|-----------|------------|----------------|
| correlation | bivariate (x, y) | m/2 · log(1 − r²) | 1 (one ρ) |
| mean | r_t | m/2 · log(ŝ²_z) | 2 (μ, σ²) |
| variance | r_t² | m/2 · log(ŝ²_z) | 2 |
| pairwise covariance | r_{i,t} · r_{j,t} | m/2 · log(ŝ²_z) | 2 |
| joint covariance matrix | R_t ∈ ℝᵖ | m/2 · log\|Σ̂\| | p(p+3)/2 + 1 |
| joint correlation matrix | R_t ∈ ℝᵖ (standardised) | m/2 · log\|Ρ̂\| | p(p+1)/2 + 1 |

where ŝ²_z = mean(z²) − mean(z)² is the biased sample variance of z over the segment, Σ̂ is the biased sample covariance matrix, and Ρ̂ is the sample correlation matrix of standardised returns.  The joint covariance and correlation matrix costs use a Cholesky-based log-determinant.

The BIC penalty is `k_params × log(n)` where `k_params = params_per_seg × m − 1` (one parameter per segment times m segments, minus one for the first segment which has no preceding changepoint location).

---

## Programs

### Correlation changepoints

**`xreturns_corr.f90`**  
Computes return statistics and the full correlation matrix. No changepoint detection; useful as a baseline.

**`xcorr_sim.f90`**  
Simulates a bivariate normal series with known piecewise-constant correlation (segments of different ρ), then detects changepoints. Prints true vs estimated parameters and writes simulated data to `sim_data.txt`.

**`xreturns_corr_cp.f90`** · **`xreturns_corr_cp.py`** · **`xreturns_corr_cp.cpp`**  
Reads `spy_efa_eem_tlt.csv`, computes percent returns, and for each asset pair finds changepoints in correlation.  After processing all pairs, prints a summary matrix:

```
             SPY   EFA   EEM   TLT
     SPY       0     2     3     1
     EFA       4     0     2     2
     EEM       5     3     0     1
     TLT       3     4     3     0
```

Above the diagonal: BIC changepoint count.  Below: AIC count.  Diagonal: 0.

The Python version (`xreturns_corr_cp.py`) uses NumPy prefix-sum vectorisation and matches the Fortran output exactly; it requires no Numba.  The C++ version uses an end-major cost matrix layout for cache efficiency.

Key parameters (compile-time constants at the top of each file):

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `MAX_CP` / `max_cp` | 20 | Maximum changepoints per pair |
| `MIN_SEG_LEN` / `min_seg_len` | 50 | Minimum observations per segment |
| `SYM_ALLOWED` / `sym_allowed` | (empty) | Subset of tickers; empty = all |
| `max_days` | 0 | Keep latest/earliest N returns; 0 = all |
| `latest` | `.true.` | If `max_days > 0`: `.true.` = most recent |

**`xcorr_fit.f90`**  
Fits the correlation changepoint model to a single pair with detailed diagnostic output.

---

### Variance and mean changepoints (univariate, per asset)

**`xreturns_variance_cp.f90`**  
Finds changepoints in the variance of each asset independently.  The derived series is z_t = r_t², and the cost is the profile normal mean-shift log-likelihood m/2 · log(ŝ²_z).  Prints a summary table of BIC/AIC changepoint counts per asset.

**`xreturns_mean_cp.f90`**  
Finds changepoints in the return mean of each asset.  Same cost function, with z_t = r_t.

**`xreturns_variance_pwl.f90`**  
Fits a **continuous piecewise-linear** path to squared returns by minimising the sum of squared residuals.  The fitted variance path v_t is continuous at the knots (changepoints) and linear between them.  This is an exact optimizer for the PWL-SSE objective, not a likelihood model.

**`xreturns_variance_pwl_fast.f90`**  
Accelerated version of `xreturns_variance_pwl.f90`.

---

### Pairwise covariance / variance changepoints

**`xreturns_cov_cp.f90`**  
For every pair (i, j) with i ≤ j (including diagonal i = j for variance), forms z_t = r_{i,t} · r_{j,t} and finds mean-shift changepoints.  The diagonal pair (i, i) detects variance changepoints; off-diagonal pairs detect covariance changepoints.  Prints a summary matrix analogous to `xreturns_corr_cp.f90`.

Note: variance of financial returns is genuinely time-varying (GARCH-like), so many changepoints are typically detected.  The normal mean-shift model is also misspecified for z_t = r² (which is non-negative and heavy-tailed).  Results are more interpretable for covariance pairs than for variance.

---

### Joint covariance matrix changepoints

**`xreturns_covmat_cp.f90`**  
Finds a **single set** of changepoints for all assets simultaneously by fitting a multivariate normal model with piecewise-constant p × p covariance matrix Σ.  The cost for a segment is m/2 · log|Σ̂|, computed via a Cholesky log-determinant.  After model selection, prints for each BIC-chosen segment the correlation lower triangle, annualised standard deviations, and optionally PCA loadings:

```
Segment 2: 2008-09-02 to 2012-06-29 (984 obs)
  Correlation:
             SPY     EFA     EEM     TLT
     SPY   1.000
     EFA   0.964   1.000
     EEM   0.932   0.956   1.000
     TLT  -0.418  -0.370  -0.310   1.000
   *SD*   0.292   0.270   0.308   0.148
```

`*SD*` is the annualised standard deviation (× √252) of unscaled returns.

Prefix sums of outer products (`SR`, `SP`) reduce cost matrix construction to O(n · p²) setup + O(n² · p²) inner work, with an O(p³) Cholesky per cell.

**`xcovar_sim.f90`**  
Simulation counterpart of `xreturns_covmat_cp.f90`.  Generates a p × n multivariate normal series with known piecewise-constant covariance matrices (specified as correlations + standard deviations per segment), runs the joint covariance changepoint algorithm, and prints true vs estimated parameters side by side.

Default setup: p = 3 assets, n = 2000 observations, 3 segments (changepoints at observations 500 and 1200) with substantially different correlation structure and volatility between segments.

---

### Joint correlation matrix changepoints

**`xreturns_corrmat_cp.f90`**  
Finds changepoints in the full p × p **correlation** matrix jointly.  Returns are first standardised so that the cost captures only changes in correlation structure, not changes in volatility.  Two normalisation schemes are selectable via `use_ewma`:

| Scheme | Description |
|--------|-------------|
| `use_ewma = .false.` | Global: subtract full-sample mean, divide by full-sample std dev.  Simple but has look-ahead bias. |
| `use_ewma = .true.` | EWMA (RiskMetrics): conditional σ_t updated at each step with decay λ = 0.94.  Removes volatility clustering before the correlation test. |

BIC parameter count per segment: p(p+1)/2 + 1 (correlations + means + changepoint location).

Additional output options (compile-time parameters):

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `print_segs` | 0 | 0 = BIC model only; 1 = 0 through BIC; 2 = all models |
| `print_diffs` | `.true.` | Bonferroni-adjusted Fisher z-tests for correlation changes between segments |
| `alpha_diff` | 0.05 | Significance level for the Fisher z-test |
| `print_pca_cov` | `.false.` | PCA of segment covariance matrix |
| `print_pca_corr` | `.true.` | PCA of segment correlation matrix |

**`xreturns_corrmat_ol.f90`**  
**Online (expanding-window)** version of `xreturns_corrmat_cp.f90`.  At each step t the DP is re-run on data 1:t only, so no future data are used.  Prints a dated table showing how the BIC-optimal changepoint dates evolve as new observations arrive:

```
as-of          n_cp  changepoints
----------------------------------------------------------------------
2009-03-31        1  2008-09-15
2009-06-30        1  2008-09-15
2009-09-30        2  2007-07-26 2008-09-15
...
```

The DP is re-run every `step` observations (default 63, ≈ 1 quarter).  Total work scales as O(n³/step).

**`xreturns_corrsub_cp.f90`**  
Runs up to four changepoint algorithms across all distinct asset subsets of specified sizes, controlled by compile-time toggles:

| Toggle | Algorithm |
|--------|-----------|
| `do_mean` | Mean-shift changepoints per asset (univariate) |
| `do_variance` | Variance changepoints per asset (on squared returns) |
| `do_cov` | Joint covariance-matrix changepoints for each subset |
| `do_corr` | Joint correlation-matrix changepoints for each subset (EWMA-standardised) |

`sub_sizes` specifies the subset sizes to analyse (default `[2, 5]`); `max_subsets` caps the number of subsets per size when C(n_assets, k) is large.  A summary of AIC/BIC changepoint counts is printed for mean and variance.

---

### Principal component analysis

**`xreturns_pca.f90`**  
Reads asset prices, computes percent returns, and performs PCA on the full-sample covariance matrix via the Jacobi eigenvalue method.  Prints the labeled covariance matrix, PC loadings, per-component variance explained, and cumulative variance explained.

---

## Modules

| File | Purpose |
|------|---------|
| `changepoint.f90` | Cost matrices (`cost_matrix`, `mean_shift_cost_matrix`, `multivar_cost_matrix`), DP solver (`solve_changepoints`), and backtracking helper (`segment_ends`) |
| `io_utils.f90` | `print_model_selection`; segment printing for corrmat, covmat, and univariate models; `print_pca_loadings`; `keep_obs` |
| `pca_jacobi.f90` | `jacobi_eigen_sym` (Jacobi eigenvalue method) and `principal_components_cov` |
| `basic_stats.f90` | `mean`, `cov_mat`, `col_stats_ignore_nan`, `standardize_returns`, `biased_cov_sd`, `print_corr_mat`, `print_acf` |
| `util.f90` | `sort_int`, `set_segment_values`, `print_wall_time`, `next_combination`, `n_choose_k`, `cumul_sum`, `print_square_matrix` |
| `sim_changepoint.f90` | `generate_series` (bivariate), `generate_multivar_series` (multivariate) |
| `dataframe_index_date.f90` | Date-indexed DataFrame type: `read_csv`, `pct_change`, `keep_rows`, `select`, etc. |
| `df_index_date_ops_mod.f90` | Index operations for the DataFrame |
| `date.f90` | `Date` type and arithmetic |
| `random.f90` | `rnorm` (scalar and vector) |
| `kind.f90` | `dp = real64`, `long_int = int64` |
| `constants.f90` | Physical and mathematical constants |

---

## Building

GNU make is required.  The compiler is gfortran (set `FC` in the Makefile to change).

```bash
# Build all programs
make -f Makefile.xcorr

# Build and run individual programs
make -f Makefile.xcorr run_sim          # xcorr_sim
make -f Makefile.xcorr run_covar_sim    # xcovar_sim
make -f Makefile.xcorr run_fit          # xcorr_fit
make -f Makefile.xcorr run_corr         # xreturns_corr_cp
make -f Makefile.xcorr run_cov          # xreturns_cov_cp
make -f Makefile.xcorr run_covmat       # xreturns_covmat_cp
make -f Makefile.xcorr run_corrmat      # xreturns_corrmat_cp
make -f Makefile.xcorr run_corrmat_ol   # xreturns_corrmat_ol
make -f Makefile.xcorr run_corrsub      # xreturns_corrsub_cp
make -f Makefile.xcorr run_pca          # xreturns_pca
make -f Makefile.xcorr run_var          # xreturns_variance_cp
make -f Makefile.xcorr run_mean         # xreturns_mean_cp
make -f Makefile.xcorr run_var_pwl      # xreturns_variance_pwl
make -f Makefile.xcorr run_var_pwl_fast # xreturns_variance_pwl_fast

# Clean
make -f Makefile.xcorr clean
```

The Python program requires NumPy and pandas and runs without compilation:

```bash
python xreturns_corr_cp.py
```

The C++ program requires a C++17 compiler:

```bash
g++ -O2 -std=c++17 -o xreturns_corr_cp xreturns_corr_cp.cpp
./xreturns_corr_cp
```

---

## Runtime

Timings on a typical desktop for the full `spy_efa_eem_tlt.csv` dataset (n ≈ 5 500 daily returns, 4 assets, `max_cp = 20`):

| Program | Time |
|---------|------|
| `xreturns_corr_cp` (Fortran) | ~2 s |
| `xreturns_corr_cp.py` (Python/NumPy) | ~10 s |
| `xreturns_corr_cp` (C++) | ~1 s |
| `xreturns_covmat_cp` | ~15 s |
| `xreturns_corrmat_cp` | ~15 s |
| `xcovar_sim` (n = 2 000) | ~0.3 s |
| `xreturns_corrmat_ol` (step = 63) | ~minutes |

The dominant cost in all programs is the O(n² × max_m) dynamic programming step.  The joint covariance and correlation matrix programs add an O(p³) Cholesky per cost matrix cell.  The online program repeats the full DP at each step, so total work is O(n³/step).

---

## References

The piecewise-constant correlation changepoint model is based on:

- Galeano, P. and Wied, D. (2014). [Multiple break detection in the correlation structure of random variables](https://arxiv.org/abs/1206.5367). *Computational Statistics & Data Analysis*, 76, 262–282.

BIC for changepoint models:

- Yao, Y.-C. (1988). [Estimating the number of change-points via Schwarz criterion](https://www.sciencedirect.com/science/article/abs/pii/0167715288901186). *Statistics & Probability Letters*, 6(3), 181–189.

Continuous piecewise-linear fitting:

- Bellman, R. and Roth, R. (1969). [Curve fitting by segmented straight lines](https://www.tandfondle.com/doi/pdf/10.1080/01621459.1969.10501038). *Journal of the American Statistical Association*, 64(327), 1079–1084.
