!> Program: xreturns_corr_cp
!!
!! Reads asset prices, computes returns, and estimates changepoints in 
!! correlation for each pair of assets using dynamic programming.
!!
!! It outputs model selection statistics (AIC, BIC) and 
!! the detected changepoint locations.

program xreturns_corr_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use basic_stats_mod, only: mean
    use changepoint_mod, only: cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time, sort_int
    use io_utils_mod, only: print_model_selection, print_estimated_parameters_dates
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp+1, min_seg_len = 50, max_assets = 1000
    integer :: n, n_col, i, j, ba, bb
    real(kind=dp), allocatable :: dp_table(:,:), cost(:,:)
    integer, allocatable :: parent(:,:), cp_aic(:,:), cp_bic(:,:)
    integer(kind=long_int) :: t_start
    character (len=10), allocatable :: date_labels(:)
    character (len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv", &
       sym_allowed(*) = [character (len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    logical, parameter :: resample_ret = .false., print_param=.true.
    integer, parameter :: max_days = 0     ! 0 = use all returns
    logical, parameter :: latest  = .true. ! if max_days>0: .true.=latest, .false.=earliest
    call system_clock(t_start)
    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read",nrow(df_px),"days and",ncol(df_px),"columns from " // trim(prices_file)
    ! Subset prices before pct_change so that n and values(2:,:) need no special handling.
    ! max_days+1 prices yield exactly max_days returns.
    if (max_days > 0) then
        df_px = df_px%keep_rows(max_days + 1, latest=latest)
        print "('using ',a,' ',i0,' returns (',a,' to ',a,')')", &
            merge('latest  ','earliest', latest), nrow(df_px) - 1, &
            trim(df_px%index(2)%to_str()), trim(df_px%index(nrow(df_px))%to_str())
    end if
    df_ret = scale_ret * df_px%pct_change()
    if (scale_ret /= 1.0_dp) print "('return scaling: ', f0.4)", scale_ret
    if (resample_ret) then
       df_ret = df_ret%resample()
       print*,"resampled returns!"
    end if
    n = size(df_ret%values, 1) - 1
    n_col = ncol(df_ret)
    date_labels = df_ret%index%to_str()    

    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n))
    allocate(cp_aic(n_col, n_col), cp_bic(n_col, n_col))
    cp_aic = 0
    cp_bic = 0

    ! Loop over each pair of assets
    do i = 1, n_col
        do j = i + 1, n_col
            print "(/,a)", "changes in " // trim(df_ret%columns(i)) // "-" // trim(df_ret%columns(j)) // " correlation"
            cost = cost_matrix(df_ret%values(2:, i), df_ret%values(2:, j), min_seg_len=min_seg_len)
            call solve_changepoints(max_m, cost, dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb)
            cp_aic(i, j) = ba   ! above-diagonal → AIC (stored lower-tri when printed below)
            cp_bic(i, j) = bb   ! above-diagonal → BIC
            if (print_param) call print_estimated_parameters_dates(max_m, parent, &
                    df_ret%values(2:, i), df_ret%values(2:, j), date_labels)
        end do
    end do

    ! Print summary matrix: above-diagonal = BIC changepoints, below = AIC
    print "(/,a)", "Summary: # changepoints chosen by BIC (above diagonal) and AIC (below diagonal)"
    ! Header row
    write(*, "(a6)", advance='no') ""
    do j = 1, n_col
        write(*, "(a6)", advance='no') trim(df_ret%columns(j))
    end do
    print *
    do i = 1, n_col
        write(*, "(a6)", advance='no') trim(df_ret%columns(i))
        do j = 1, n_col
            if (i == j) then
                write(*, "(i6)", advance='no') 0
            else if (j > i) then
                write(*, "(i6)", advance='no') cp_bic(i, j)
            else
                write(*, "(i6)", advance='no') cp_aic(j, i)
            end if
        end do
        print *
    end do

    deallocate(cost, dp_table, parent, cp_aic, cp_bic)
    print "(/,a)", "(3) finished xreturns_corr_cp.f90"
    call print_wall_time(t_start)
end program xreturns_corr_cp
