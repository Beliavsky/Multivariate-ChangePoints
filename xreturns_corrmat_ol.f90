!> Program: xreturns_corrmat_ol
!!
!! Online (expanding-window) changepoint detection for the full p×p correlation
!! matrix.  At each step t the DP is re-run on data 1:t only, so no future data
!! are used.  The BIC-optimal number of changepoints and their dates are printed
!! as a table showing how the segmentation evolves over time.
!!
!! Performance note: at each step the full t×t cost matrix is recomputed, so
!! total work scales as O(n³/step).  For large n, increase step accordingly.
!!
!! EWMA standardization is causal (σ²_t uses only data up to t-1) so
!! standardizing on the full sample is valid here.  Global standardization
!! (use_ewma=.false.) uses the full-sample mean and std — this has mild
!! look-ahead bias and is best avoided in online use.

program xreturns_corrmat_ol
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: multivar_cost_matrix, solve_changepoints, segment_ends
    use basic_stats_mod, only: standardize_returns
    use util_mod, only: print_wall_time
    use io_utils_mod, only: keep_obs
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20         ! max changepoints considered
    integer, parameter :: max_m  = max_cp + 1 ! max segments
    integer, parameter :: min_seg_len = 50    ! minimum observations per segment
    integer, parameter :: max_assets = 1000   ! max columns read from CSV
    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv"
    character(len=5),  parameter :: sym_allowed(*) = [character(len=5) ::]  ! empty = use all
    real(kind=dp), parameter :: scale_ret   = 100.0_dp   ! returns multiplied by this (basis: 100 = percent)
    integer,       parameter :: max_days    = 0           ! 0 = use all returns
    logical,       parameter :: latest      = .true.      ! if max_days > 0: use latest observations
    logical,       parameter :: use_ewma    = .true.      ! .true. = EWMA, .false. = global standardization
    real(kind=dp), parameter :: ewma_lambda = 0.94_dp    ! EWMA decay factor (RiskMetrics daily default)
    integer,       parameter :: step  = 63                ! re-run every step obs (~1 quarter)
    integer,       parameter :: min_n = 2 * min_seg_len  ! minimum obs before first run
    logical,       parameter :: resample_ret = .false.   ! .true. = shuffle rows (null hypothesis check)

    integer :: n, n_col, pps, best_bic, t, j
    real(kind=dp), allocatable :: dp_table_t(:,:), R(:,:), R_std(:,:)
    integer,       allocatable :: parent_t(:,:), cp_ends(:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: ret_dates(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read", nrow(df_px), "days and", ncol(df_px), &
        "columns from " // trim(prices_file)
    if (max_days > 0) call keep_obs(df_px, max_days, latest, verbose=.true.)
    df_ret = scale_ret * df_px%pct_change(dropna=.true.)
    if (resample_ret) then
        df_ret = df_ret%resample()
        print *, "resampled returns!"
    end if

    n     = size(df_ret%values, 1)
    n_col = ncol(df_ret)
    pps   = n_col * (n_col + 1) / 2 + 1  ! params per segment: p*(p+1)/2 correlations + 1 location

    ret_dates = df_ret%index%to_str()
    R         = df_ret%values
    R_std       = standardize_returns(R, use_ewma, ewma_lambda)

    print "(/,a)", "Online changepoint detection (expanding window)"
    print "('step = ',i0,' obs,  min_n = ',i0,' obs,  p = ',i0,' assets')", step, min_n, n_col
    print "('params_per_seg = ',i0,'  (p*(p+1)/2 + 1)')", pps
    print *
    write(*, "(a12, a6, 2x, a)") "as-of", "n_cp", "changepoints"
    write(*, "(a)") repeat('-', 70)

    t = min_n - step
    do
        t = min(t + step, n)

        allocate(dp_table_t(t, max_m), parent_t(t, max_m))
        call solve_changepoints(max_m, &
            multivar_cost_matrix(R_std(1:t,:), min_seg_len=min_seg_len), &
            dp_table_t, parent_t)

        best_bic = bic_cp(dp_table_t, pps)

        write(*, "(a12, i6, 2x)", advance='no') ret_dates(t), best_bic
        if (best_bic > 0) then
            cp_ends = segment_ends(parent_t, best_bic + 1)
            do j = 1, best_bic
                write(*, "(a10, 1x)", advance='no') ret_dates(cp_ends(j))
            end do
        else
            write(*, "(a)", advance='no') "(none)"
        end if
        print *

        deallocate(dp_table_t, parent_t)
        if (allocated(cp_ends)) deallocate(cp_ends)
        if (t == n) exit
    end do

    deallocate(R, R_std)
    call print_wall_time(t_start)

contains

    pure function bic_cp(dp_tab, pps_) result(best_m_minus_1)
        !> Return the BIC-optimal number of changepoints.
        real(kind=dp), intent(in) :: dp_tab(:,:)
        integer,       intent(in) :: pps_
        integer :: best_m_minus_1, m, nt, max_m_loc
        real(kind=dp) :: bic, min_bic
        nt        = size(dp_tab, 1)
        max_m_loc = size(dp_tab, 2)
        min_bic   = huge(1.0_dp)
        best_m_minus_1 = 0
        do m = 1, max_m_loc
            if (dp_tab(nt, m) >= 1.0e19_dp) cycle
            bic = 2.0_dp * dp_tab(nt, m) + real(pps_ * m - 1, dp) * log(real(nt, dp))
            if (bic < min_bic) then
                min_bic        = bic
                best_m_minus_1 = m - 1
            end if
        end do
    end function bic_cp

end program xreturns_corrmat_ol
