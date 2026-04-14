"""
xreturns_corr_cp.py

Python equivalent of xreturns_corr_cp.f90.
Reads asset prices, computes percent returns, and for each asset pair finds
changepoints in correlation using O(n^2 * max_m) dynamic programming.

No numba -- speed is kept manageable by limiting MAX_CP (default 2).
Cost matrix uses O(n) Python iterations with numpy-vectorised inner work.
"""

import time
import numpy as np
import pandas as pd
from itertools import combinations

# ── parameters ───────────────────────────────────────────────────────────────
PRICES_FILE  = "spy_efa_eem_tlt.csv"
SCALE_RET    = 100.0
MAX_CP       = 10          # maximum changepoints per pair  ← main speed knob
MIN_SEG_LEN  = 50
SYM_ALLOWED  = None # ["SPY", "TLT"]  # subset of columns; [] or None → all columns
PRINT_PARAMS = False            # print per-segment correlation for each model
# ─────────────────────────────────────────────────────────────────────────────


def build_cost_matrix(x: np.ndarray, y: np.ndarray,
                      min_seg_len: int = 50) -> np.ndarray:
    """
    n×n cost matrix where cost[i, j] = 0.5*(j-i+1)*log(1-r²) for segment
    i..j (0-indexed, inclusive); 1e20 where j-i+1 < min_seg_len or j < i.

    Uses cumulative sums (O(n) prefix-sum setup) and O(n) Python loop with
    O(n) numpy work per iteration -- no numba required.
    Memory: n*n float64 (e.g. ~200 MB for n=5000).
    """
    n    = len(x)
    cost = np.full((n, n), 1e20, dtype=np.float64)

    # length-(n+1) prefix sums; s[0]=0, s[k] = sum(x[0..k-1])
    sx  = np.empty(n + 1); sx[0]  = 0.0; np.cumsum(x,   out=sx[1:])
    sy  = np.empty(n + 1); sy[0]  = 0.0; np.cumsum(y,   out=sy[1:])
    sxx = np.empty(n + 1); sxx[0] = 0.0; np.cumsum(x*x, out=sxx[1:])
    syy = np.empty(n + 1); syy[0] = 0.0; np.cumsum(y*y, out=syy[1:])
    sxy = np.empty(n + 1); sxy[0] = 0.0; np.cumsum(x*y, out=sxy[1:])

    for i in range(n - min_seg_len + 1):
        # j1: prefix-sum end indices for all valid j >= i+min_seg_len-1
        j1    = np.arange(i + min_seg_len, n + 1)  # j1 = j+1 in 0-indexed obs
        seg_n = (j1 - i).astype(np.float64)        # segment lengths

        ma = (sx[j1]  - sx[i])  / seg_n
        mb = (sy[j1]  - sy[i])  / seg_n
        va = (sxx[j1] - sxx[i]) / seg_n - ma**2   # biased variance of x
        vb = (syy[j1] - syy[i]) / seg_n - mb**2   # biased variance of y

        r = ((sxy[j1] - sxy[i]) / seg_n - ma * mb) / \
            np.sqrt(np.maximum(va * vb, 1e-20))
        np.clip(r, -(1.0 - 1e-10), 1.0 - 1e-10, out=r)

        cost[i, i + min_seg_len - 1:] = \
            0.5 * seg_n * np.log(np.maximum(1.0 - r**2, 1e-10))

    return cost


def solve_changepoints(cost: np.ndarray,
                       max_m: int) -> tuple[np.ndarray, np.ndarray]:
    """
    dp[i, m] = minimum total cost for observations 0..i split into m+1
               segments (0-indexed m; Fortran m+1 segments).
    par[i, m] = 0-indexed last position of segment m (the predecessor boundary).

    Direct Python translation of the Fortran solve_changepoints, with numpy
    slice ops replacing the innermost k-loop.
    """
    n   = cost.shape[0]
    dp  = np.full((n, max_m), 1e20, dtype=np.float64)
    par = np.zeros((n, max_m), dtype=np.int32)

    dp[:, 0] = cost[0, :]   # m=0 (1 segment): cost from obs 0 to i

    for m in range(1, max_m):            # Fortran: m = 2..max_m
        for i in range(m, n):            # Fortran: i = m..n
            # k runs from m-1 to i-1 (0-indexed); last segment is k+1..i
            dp_prev  = dp[m - 1 : i, m - 1]  # dp[k, m-1] for k = m-1..i-1
            cost_seg = cost[m : i + 1, i]     # cost[k+1, i] for k = m-1..i-1
            vals     = dp_prev + cost_seg
            best     = int(np.argmin(vals))
            dp[i, m]  = vals[best]
            par[i, m] = (m - 1) + best        # k = (m-1) + best_idx

    return dp, par


def reconstruct_cps(par: np.ndarray, n: int, m_segs: int) -> list[int]:
    """
    Returns sorted list of 0-indexed last positions of segments 1..m_segs-1
    (i.e., the m_segs-1 changepoint locations as end-of-segment indices).
    Fortran equivalent prints these as 1-indexed: add 1 to each value.
    """
    cp  = n - 1
    cps = []
    for k in range(m_segs - 1, 0, -1):
        cp = int(par[cp, k])
        cps.append(cp)
    return sorted(cps)


def seg_corr(x: np.ndarray, y: np.ndarray, start: int, end: int) -> float:
    """Pearson correlation for x[start:end+1]."""
    sx, sy = x[start:end+1], y[start:end+1]
    dx, dy = sx - sx.mean(), sy - sy.mean()
    den = np.sqrt((dx**2).sum() * (dy**2).sum())
    return float((dx * dy).sum() / den) if den > 1e-15 else 0.0


def print_pair_results(dp: np.ndarray, par: np.ndarray,
                       x: np.ndarray, y: np.ndarray,
                       n: int, max_m: int,
                       dates: list[str],
                       print_params: bool = False) -> None:
    """Mirrors print_model_selection + print_estimated_parameters_dates."""
    print(f"{'M':>2}  {'LL':>10}  {'AIC':>10}  {'BIC':>10}  ChangePoints")

    aic_vals = np.full(max_m, np.inf)
    bic_vals = np.full(max_m, np.inf)
    log_n    = np.log(n)

    for m_segs in range(1, max_m + 1):
        m          = m_segs - 1
        total_cost = dp[n - 1, m]
        if total_cost >= 1e19:
            continue
        ll       = -total_cost
        k_params = 2 * m_segs - 1
        aic      = -2.0 * ll + 2.0 * k_params
        bic      = -2.0 * ll + k_params * log_n
        aic_vals[m] = aic
        bic_vals[m] = bic

        cps = reconstruct_cps(par, n, m_segs)
        # Print 1-indexed positions matching Fortran, then dates in parentheses
        cp_str = "  ".join(
            f"{c + 1} ({dates[c]})" for c in cps
        )
        print(f"{m_segs:>2}  {ll:>10.2f}  {aic:>10.2f}  {bic:>10.2f}  {cp_str}")

        if print_params:
            bounds = [-1] + cps + [n - 1]
            print(f"  {'Start':>12}  {'End':>12}  {'#obs':>5}  {'Corr':>8}")
            for k in range(m_segs):
                s = bounds[k] + 1
                e = bounds[k + 1]
                print(f"  {dates[s]:>12}  {dates[e]:>12}  "
                      f"{e - s + 1:>5}  {seg_corr(x, y, s, e):>8.4f}")

    best_aic = int(np.argmin(aic_vals)) + 1
    best_bic = int(np.argmin(bic_vals)) + 1
    print(f"Changepoints chosen by AIC, BIC: {best_aic - 1}  {best_bic - 1}")


def main() -> None:
    t0 = time.time()

    # ── load prices and compute returns ──────────────────────────────────────
    df_px = pd.read_csv(PRICES_FILE, index_col=0, parse_dates=True)
    if SYM_ALLOWED:
        df_px = df_px[[c for c in SYM_ALLOWED if c in df_px.columns]]
    cols = list(df_px.columns)

    df_ret = SCALE_RET * df_px.pct_change().iloc[1:]
    n      = len(df_ret)
    dates  = [d.strftime("%Y-%m-%d") for d in df_ret.index]
    ret    = df_ret.to_numpy(dtype=np.float64)

    print(f"read {n + 1} days and {len(cols)} columns from {PRICES_FILE}")
    if SCALE_RET != 1.0:
        print(f"return scaling: {SCALE_RET:.4f}")

    max_m = MAX_CP + 1

    # ── loop over asset pairs ────────────────────────────────────────────────
    for i, j in combinations(range(len(cols)), 2):
        pair = f"{cols[i]}-{cols[j]}"
        print(f"\nchanges in {pair} correlation")
        t_pair = time.time()
        cost    = build_cost_matrix(ret[:, i], ret[:, j], MIN_SEG_LEN)
        dp, par = solve_changepoints(cost, max_m)
        print_pair_results(dp, par, ret[:, i], ret[:, j],
                           n, max_m, dates, PRINT_PARAMS)
        print(f"  ({time.time() - t_pair:.1f}s)")

    print(f"\nwall time elapsed: {time.time() - t0:.2f}s")


if __name__ == "__main__":
    main()
