!> Program: xreturns_corrmat_cp
!!
!! Reads asset prices, computes returns, and finds changepoints in the FULL
!! p×p correlation matrix jointly using dynamic programming.
!!
!! Model: within each segment the standardised returns are i.i.d. multivariate
!! normal with correlation matrix Rho.  Returns are first standardised so that
!! the multivariate cost function captures only changes in correlation structure,
!! not changes in volatility.  Two normalization schemes are available via the
!! use_ewma parameter:
!!
!!   use_ewma = .false.  Global standardisation: subtract full-sample mean and
!!                       divide by full-sample std dev for each asset.  Simple
!!                       but has look-ahead bias and assumes constant volatility.
!!
!!   use_ewma = .true.   EWMA (RiskMetrics) standardisation: at each time step t
!!                       the conditional std dev is
!!                           sigma_t = sqrt(lambda*sigma_{t-1}^2 + (1-lambda)*r_{t-1}^2)
!!                       initialised with the full-sample variance.  This removes
!!                       volatility clustering before the correlation test, giving
!!                       a cleaner separation of correlation changes from variance
!!                       changes.  lambda = 0.94 is the RiskMetrics daily default.
!!
!! The profile log-likelihood cost for a segment is
!!
!!     cost(i..j) = m/2 * log|Rho_hat(i:j)|
!!
!! computed via multivar_cost_matrix on the standardised returns.
!!
!! BIC parameter count per segment: p*(p-1)/2 off-diagonal correlations
!! + p means + 1 changepoint location = p*(p+1)/2 + 1.
!!
!! After model selection the estimated correlation matrix is printed for each
!! segment of the BIC-chosen model, together with annualised standard deviations
!! computed from the original (unstandardised) returns in that segment.

program xreturns_corrmat_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: multivar_cost_matrix, solve_changepoints, segment_ends
    use basic_stats_mod, only: standardize_returns
    use util_mod, only: print_wall_time
    use io_utils_mod, only: print_model_selection, print_corrmat_model, corrmat_model_range, &
        keep_obs
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20          ! max changepoints considered
    integer, parameter :: max_m = max_cp + 1   ! max segments = max_cp + 1
    integer, parameter :: min_seg_len = 50     ! minimum observations per segment
    integer, parameter :: max_assets = 1000    ! max columns read from CSV
    character(len=*), parameter :: prices_file = "asset_class_etf_prices.csv" ! "spy_efa_eem_tlt.csv"
    character(len=5),  parameter :: sym_allowed(*) = [character(len=5) :: "SPY", "USO"]  ! empty = use all
    real(kind=dp), parameter :: scale_ret = 100.0_dp  ! returns multiplied by this (basis: 100 = percent)
    integer, parameter :: max_days = 0     ! 0 = use all returns
    logical, parameter :: latest  = .true. ! if max_days>0: .true.=latest, .false.=earliest
    logical, parameter :: use_ewma = .true.          ! .true. = EWMA, .false. = global standardization
    logical, parameter :: print_all_models = .false. ! print AIC/BIC table for all numbers of changepoints
    real(kind=dp), parameter :: ewma_lambda = 0.94_dp ! EWMA decay factor (RiskMetrics daily default)
    ! print_segs: 0 = BIC-chosen model only
    !             1 = 0 changepoints through BIC-chosen
    !             2 = all models studied
    integer, parameter :: print_segs = 0
    logical, parameter :: print_diffs = .true.   ! test sig. corr. changes between consecutive segments
    real(kind=dp), parameter :: alpha_diff = 0.05_dp ! significance level for the Fisher z-test that compares correlations between consecutive segments
    logical, parameter :: print_pca_cov  = .false. ! PCA of segment covariance matrix
    logical, parameter :: print_pca_corr = .true.  ! PCA of segment correlation matrix

    integer :: n, n_col, pps, best_aic, best_bic, ms, m_lo, m_hi
    real(kind=dp), allocatable :: dp_table(:,:), R(:,:), R_std(:,:)
    integer, allocatable :: parent(:,:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: date_labels(:), ret_dates(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read", nrow(df_px), "days and", ncol(df_px), &
        "columns from " // trim(prices_file)
    if (max_days > 0) call keep_obs(df_px, max_days, latest, verbose=.true.)
    df_ret = scale_ret * df_px%pct_change()
    if (scale_ret /= 1.0_dp) print "('return scaling: ', f0.4)", scale_ret

    n     = size(df_ret%values, 1) - 1
    n_col = ncol(df_ret)

    ! ret_dates(i) is the date of return i (1-indexed), skipping the NaN first row
    date_labels = df_ret%index%to_str()
    ret_dates = date_labels(2:n+1)

    ! R(i, a) = return of asset a on day i (original, unscaled)
    R = df_ret%values(2:, :)

    R_std = standardize_returns(R, use_ewma, ewma_lambda)

    ! ── cost matrix and DP on standardised returns ───────────────────────────────
    print "(/,'building ',i0,'x',i0,' joint cost matrix for ',i0,' assets...')", &
        n, n, n_col
    allocate(dp_table(n, max_m), parent(n, max_m))

    print "(a)", "running DP..."
    call solve_changepoints(max_m, multivar_cost_matrix(R_std, min_seg_len=min_seg_len), &
         dp_table, parent)
    deallocate(R_std)

    ! ── model selection ──────────────────────────────────────────────────────────
    ! k = pps*m - 1,  pps = p*(p-1)/2 + p + 1 = p*(p+1)/2 + 1
    pps = n_col * (n_col + 1) / 2 + 1
    print "('params_per_seg = ',i0,'  (p=',i0,', p*(p+1)/2+1 = ',i0,')')", pps, n_col, pps
    call print_model_selection(dp_table, parent, best_aic_cp=best_aic, &
       best_bic_cp=best_bic, params_per_seg=pps, print_each=print_all_models)
    call corrmat_model_range(print_segs, best_bic, max_m, m_lo, m_hi) ! set m_lo and m_hi

    do ms = m_lo, m_hi
        if (dp_table(n, ms) >= 1.0e19_dp) cycle
        call print_corrmat_model(best_bic, segment_ends(parent, ms), R, df_ret%columns, &
           ret_dates, scale_ret, print_diffs, alpha_diff, &
           do_pca_cov=print_pca_cov, do_pca_corr=print_pca_corr)
    end do

    deallocate(dp_table, parent, R)
    call print_wall_time(t_start)

end program xreturns_corrmat_cp
