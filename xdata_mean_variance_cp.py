"""
xdata_mean_variance_cp.py

Reads a matrix data file and for each column finds changepoints in the mean
and/or variance using a profile Gaussian (log-variance) cost function.

Input file format:
  First line:  # nrow ncol [any other text]
  Further lines beginning with '#' are treated as comments and skipped.
  Remaining lines: nrow rows of ncol whitespace-separated real values.

Output matches xdata_mean_variance_cp.f90.
"""

import sys
import time
import numpy as np

try:
    import ruptures as rpt
    from ruptures.base import BaseCost
except ImportError:
    sys.exit("ruptures not installed:  pip install ruptures")

# ── parameters ────────────────────────────────────────────────────────────────
DATA_FILE   = "sim_matrix_300_obs.txt"
DO_MEAN     = False
DO_VARIANCE = True
USE_LOG_VAR = True
VAR_OFFSET  = 0.01      # c in log(c + x^2); ignored when USE_LOG_VAR=False
MIN_SEG_LEN = 20
MAX_COL     = 3         # 0 = use all columns
MAX_CP      = 100
PRINT_TABLE = False     # True  => print M/cost/AIC/BIC table for each series
PRINT_SEGS  = 1         # 0=none, 1=BIC model segments
PPS         = 3         # parameters per segment (mu + sigma^2 + 1 location/break)
# ──────────────────────────────────────────────────────────────────────────────


def read_matrix(filename, ncol_max=0):
    """Read matrix file with '# nrow ncol ...' header line."""
    nrow = ncol = 0
    data_lines = []
    with open(filename) as f:
        first_comment = True
        for line in f:
            if line.startswith('#'):
                if first_comment:
                    parts = line.lstrip('#').split()
                    nrow, ncol = int(parts[0]), int(parts[1])
                    first_comment = False
            else:
                data_lines.append(line)
    nc = min(ncol, ncol_max) if ncol_max > 0 else ncol
    x = np.empty((nrow, nc))
    for i, line in enumerate(data_lines[:nrow]):
        x[i] = list(map(float, line.split()[:nc]))
    print(f"read {nrow} rows x {nc} cols from {filename}")
    return x


class GaussianProfileCost(BaseCost):
    """
    Profile Gaussian cost per segment: n/2 * log(sample_variance).
    Matches the Fortran mean_shift_cost_matrix used in solve_changepoints.
    """
    model = "gaussian_profile"
    min_size = 2

    def fit(self, signal):
        self.signal = signal.ravel()
        return self

    def error(self, start, end):
        seg = self.signal[start:end]
        m = len(seg)
        if m < 2:
            return 1e20
        v = np.var(seg, ddof=1)
        return 1e20 if v <= 0 else 0.5 * m * np.log(v)


def total_cost(z, breakpoints):
    """Sum of Gaussian profile costs over all segments (breakpoints ends with n)."""
    total, prev = 0.0, 0
    for bp in breakpoints:
        seg = z[prev:bp]
        m = len(seg)
        if m >= 2:
            v = np.var(seg, ddof=1)
            if v > 0:
                total += 0.5 * m * np.log(v)
        prev = bp
    return total


def bic(cost, n_segs, n):
    return 2.0 * cost + (PPS * n_segs - 1) * np.log(n)


def aic(cost, n_segs):
    return 2.0 * cost + 2.0 * (PPS * n_segs - 1)


def print_segments(x, breakpoints, label):
    """Print segment statistics: mean, sd, min, max, first, last."""
    n_bkps = len(breakpoints) - 1
    print(f"\nBIC-selected model ({n_bkps} changepoint(s)) -- {label}")
    print(f"{'start':>10s}{'end':>10s}{'n':>8s}"
          f"{'mean':>12s}{'sd':>12s}{'min':>12s}{'max':>12s}{'first':>12s}{'last':>12s}")
    prev = 0
    for bp in breakpoints:
        seg = x[prev:bp]
        m = len(seg)
        mean = np.mean(seg)
        sd   = np.std(seg, ddof=1) if m > 1 else 0.0
        print(f"{prev+1:>10d}{bp:>10d}{m:>8d}"
              f"{mean:>12.4f}{sd:>12.4f}{np.min(seg):>12.4f}"
              f"{np.max(seg):>12.4f}{seg[0]:>12.4f}{seg[-1]:>12.4f}")
        prev = bp

def analyze(z, label, x_orig=None, jump=1):
    """
    Detect changepoints in z.  Returns (bic_n_bkps, aic_n_bkps, bkps_for_bic).
    When PRINT_SEGS > 0, segment stats are printed using x_orig (defaults to z).
    """
    n = len(z)
    max_bkps = min(MAX_CP, n // MIN_SEG_LEN - 1)
    if x_orig is None:
        x_orig = z

    # Dynp finds the optimal segmentation for each fixed number of breakpoints.
    # jump=1 considers every observation as a candidate breakpoint (exact, like Fortran DP).
    algo = rpt.Dynp(custom_cost=GaussianProfileCost(), min_size=MIN_SEG_LEN, jump=jump)
    algo.fit(z.reshape(-1, 1))

    if PRINT_TABLE:
        print(f"\n{'m':>4s}{'cost':>12s}{'AIC':>12s}{'BIC':>12s}    ChangePoints")

    best_bic_val = best_aic_val = np.inf
    bkps_bic = bkps_aic = [n]
    bic_nb = aic_nb = 0

    for k in range(0, max_bkps + 1):
        try:
            bkps = [n] if k == 0 else algo.predict(n_bkps=k)
        except Exception:
            break
        cost  = total_cost(z, bkps)
        bic_v = bic(cost, k + 1, n)
        aic_v = aic(cost, k + 1)

        if PRINT_TABLE:
            cp_str = " ".join(str(b) for b in bkps[:-1])
            print(f"{k+1:>4d}{cost:>12.2f}{aic_v:>12.2f}{bic_v:>12.2f}    {cp_str}")

        if bic_v < best_bic_val:
            best_bic_val = bic_v;  bic_nb = k;  bkps_bic = bkps
        if aic_v < best_aic_val:
            best_aic_val = aic_v;  aic_nb = k;  bkps_aic = bkps

    print(f"\nChangepoints chosen by AIC, BIC: {aic_nb} {bic_nb}\n")

    if PRINT_SEGS > 0:
        print_segments(x_orig, bkps_bic, label)

    return bic_nb, aic_nb


def main():
    t0 = time.time()
    x = read_matrix(DATA_FILE, ncol_max=MAX_COL)
    n, n_col = x.shape
    col_labels = [f"x{j+1}" for j in range(n_col)]
    jump = 1
    print("jump:", jump)
    mean_bic_arr = np.zeros(n_col, dtype=int)
    mean_aic_arr = np.zeros(n_col, dtype=int)
    var_bic_arr  = np.zeros(n_col, dtype=int)
    var_aic_arr  = np.zeros(n_col, dtype=int)

    if DO_MEAN:
        for j, label in enumerate(col_labels):
            print(f"\nchanges in mean of {label}")
            bb, ba = analyze(x[:, j], label, jump=jump)
            mean_bic_arr[j], mean_aic_arr[j] = bb, ba

    if DO_VARIANCE:
        var_name = "log-variance" if USE_LOG_VAR else "variance"
        for j, label in enumerate(col_labels):
            z = np.log(VAR_OFFSET + x[:, j]**2) if USE_LOG_VAR else x[:, j]**2
            print(f"\nchanges in {var_name} of {label}")
            bb, ba = analyze(z, label, x_orig=x[:, j], jump=jump)
            var_bic_arr[j], var_aic_arr[j] = bb, ba

    # ── summary ───────────────────────────────────────────────────────────────
    if DO_MEAN:
        print("\nSummary: # mean changepoints chosen by AIC and BIC")
        print(f"{'Series':>10s}{'AIC':>8s}{'BIC':>8s}")
        for j, label in enumerate(col_labels):
            print(f"{label:>10s}{mean_aic_arr[j]:>8d}{mean_bic_arr[j]:>8d}")

    if DO_VARIANCE:
        var_name = "log-variance" if USE_LOG_VAR else "variance"
        print(f"\nSummary: # {var_name} changepoints chosen by AIC and BIC")
        print(f"{'Series':>10s}{'AIC':>8s}{'BIC':>8s}")
        for j, label in enumerate(col_labels):
            print(f"{label:>10s}{var_aic_arr[j]:>8d}{var_bic_arr[j]:>8d}")

    print(f"\nwall time elapsed (s): {time.time()-t0:.4f}")


if __name__ == "__main__":
    main()
