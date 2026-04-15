!> Program: xreturns_variance_cp
!!
!! Reads asset prices, computes returns, and for each asset finds changepoints
!! in the variance using a profile normal mean-shift model applied to the
!! squared-return series z_t = r_t^2.
!!
!! Cost function per segment: m/2 * log(sample_variance(z))
!! This profiles over both the mean and variance of z within each segment,
!! so it detects any shift in the variance process regardless of cause.
!!
!! Outputs: for each asset, a table of M, LL, AIC, BIC, changepoint positions.
!! If print_segs > 0, also prints date ranges and statistics for each segment
!! of the BIC-chosen model (and optionally a range of models via print_segs).
!! At the end: a summary table of # changepoints chosen by AIC and BIC.
!!
!! print_segs: 0 = AIC/BIC table only
!!             1 = segment details for BIC-chosen model only
!!             2 = segment details for 0 changepoints through BIC-chosen
!!             3 = segment details for all models studied

program xreturns_variance_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints, segment_ends
    use util_mod, only: print_wall_time
    use io_utils_mod, only: print_model_selection, print_return_segments, corrmat_model_range, keep_obs
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp + 1, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = "asset_class_etf_prices.csv", &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    integer, parameter :: max_days = 0       ! 0 = use all returns
    logical, parameter :: latest  = .true.   ! if max_days>0: .true.=latest, .false.=earliest
    integer, parameter :: print_segs = 1     ! 0=none, 1=BIC only, 2=0..BIC, 3=all

    integer :: n, n_col, j, ba, bb, ms, m_lo, m_hi
    real(kind=dp), allocatable :: dp_table(:,:), R(:,:), z(:)
    integer, allocatable :: parent(:,:), cp_aic(:), cp_bic(:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: ret_dates(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read", nrow(df_px), "days and", ncol(df_px), &
        "columns from " // trim(prices_file)
    if (max_days > 0) call keep_obs(df_px, max_days, latest, verbose=.true.)
    df_ret = scale_ret * df_px%pct_change(dropna=.true.)
    if (scale_ret /= 1.0_dp) print "('return scaling: ', f0.4)", scale_ret

    n     = size(df_ret%values, 1)
    n_col = ncol(df_ret)

    ret_dates = df_ret%index%to_str()

    allocate(R(n, n_col), dp_table(n, max_m), parent(n, max_m), z(n))
    R = df_ret%values
    allocate(cp_aic(n_col), cp_bic(n_col))
    cp_aic = 0
    cp_bic = 0

    do j = 1, n_col
        z = df_ret%values(:, j)**2

        print "(/,a)", "changes in variance of " // trim(df_ret%columns(j))

        call solve_changepoints(max_m, mean_shift_cost_matrix(z, min_seg_len=min_seg_len), &
            dp_table, parent)
        call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, &
            params_per_seg=3)
        cp_aic(j) = ba
        cp_bic(j) = bb

        if (print_segs > 0) then
            call corrmat_model_range(print_segs - 1, bb, max_m, m_lo, m_hi)
            do ms = m_lo, m_hi
                if (dp_table(n, ms) >= 1.0e19_dp) cycle
                call print_return_segments(bb, segment_ends(parent, ms), R(:,j), &
                    trim(df_ret%columns(j)), ret_dates, scale_ret)
            end do
        end if
    end do

    print "(/,a)", "Summary: # variance changepoints chosen by AIC and BIC"
    print "(a10,2a8)", "Series", "AIC", "BIC"
    do j = 1, n_col
        print "(a10,2i8)", trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
    end do

    deallocate(R, dp_table, parent, z, cp_aic, cp_bic)
    call print_wall_time(t_start)
end program xreturns_variance_cp
