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
!! At the end: a summary table of # changepoints chosen by AIC and BIC.

program xreturns_variance_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time
    use io_utils_mod, only: print_model_selection
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp + 1, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv", &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    integer, parameter :: max_days = 260     ! 0 = use all returns
    logical, parameter :: latest  = .true. ! if max_days>0: .true.=latest, .false.=earliest

    integer :: n, n_col, j, ba, bb
    real(kind=dp), allocatable :: dp_table(:,:), cost(:,:), z(:)
    integer, allocatable :: parent(:,:), cp_aic(:), cp_bic(:)
    integer(kind=long_int) :: t_start

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read", nrow(df_px), "days and", ncol(df_px), &
        "columns from " // trim(prices_file)
    if (max_days > 0) then
        df_px = df_px%keep_rows(max_days + 1, latest=latest)
        print "('using ',a,' ',i0,' returns (',a,' to ',a,')')", &
            merge('latest  ','earliest', latest), nrow(df_px) - 1, &
            trim(df_px%index(2)%to_str()), trim(df_px%index(nrow(df_px))%to_str())
    end if
    df_ret = scale_ret * df_px%pct_change()
    if (scale_ret /= 1.0_dp) print "('return scaling: ', f0.4)", scale_ret

    n     = size(df_ret%values, 1) - 1
    n_col = ncol(df_ret)

    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n), z(n))
    allocate(cp_aic(n_col), cp_bic(n_col))
    cp_aic = 0
    cp_bic = 0

    do j = 1, n_col
        z = df_ret%values(2:, j)**2

        print "(/,a)", "changes in variance of " // trim(df_ret%columns(j))

        cost = mean_shift_cost_matrix(z, min_seg_len=min_seg_len)
        call solve_changepoints(max_m, cost, dp_table, parent)
        call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, params_per_seg=3)
        cp_aic(j) = ba
        cp_bic(j) = bb
    end do

    print "(/,a)", "Summary: # variance changepoints chosen by AIC and BIC"
    print "(a10,2a8)", "Series", "AIC", "BIC"
    do j = 1, n_col
        print "(a10,2i8)", trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
    end do

    deallocate(cost, dp_table, parent, z, cp_aic, cp_bic)
    call print_wall_time(t_start)
end program xreturns_variance_cp
