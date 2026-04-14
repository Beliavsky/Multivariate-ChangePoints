!> Program: xreturns_cov_cp
!!
!! Reads asset prices, computes returns, and for each pair of assets (i <= j)
!! finds changepoints in the covariance (i /= j) or variance (i == j) using
!! a profile normal mean-shift model applied to the product series z_t = r_it * r_jt.
!!
!! Cost function per segment: m/2 * log(sample_variance(z))
!! This profiles over both the mean and variance of z within each segment,
!! so it detects any shift in covariance or variance regardless of cause.
!!
!! Outputs: for each pair, a table of M, LL, AIC, BIC, changepoint positions.
!! At the end: summary matrices (BIC above diagonal, AIC below) of # changepoints.

program xreturns_cov_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time
    use io_utils_mod, only: print_model_selection, print_estimated_parameters_dates
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp + 1, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv", &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    logical, parameter :: print_param = .false.
    integer, parameter :: max_days = 260     ! 0 = use all returns
    logical, parameter :: latest  = .true. ! if max_days>0: .true.=latest, .false.=earliest

    integer :: n, n_col, j1, j2, ba, bb
    real(kind=dp), allocatable :: dp_table(:,:), cost(:,:), z(:)
    integer, allocatable :: parent(:,:), cp_aic(:,:), cp_bic(:,:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: date_labels(:)
    character(len=64) :: pair_label

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
    date_labels = df_ret%index%to_str()

    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n), z(n))
    allocate(cp_aic(n_col, n_col), cp_bic(n_col, n_col))
    cp_aic = 0
    cp_bic = 0

    do j1 = 1, n_col
        do j2 = j1, n_col
            ! z_t = r_{j1,t} * r_{j2,t}
            ! diagonal (j1==j2): z_t = r²_{j1,t}  → variance changepoints
            ! off-diagonal:       z_t = r_{j1,t}*r_{j2,t} → covariance changepoints
            z = df_ret%values(2:, j1) * df_ret%values(2:, j2)

            if (j1 == j2) then
                pair_label = "variance of " // trim(df_ret%columns(j1))
            else
                pair_label = "covariance of " // trim(df_ret%columns(j1)) // &
                             "-" // trim(df_ret%columns(j2))
            end if
            print "(/,a)", "changes in " // trim(pair_label)

            cost = mean_shift_cost_matrix(z, min_seg_len=min_seg_len)
            call solve_changepoints(max_m, cost, dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, params_per_seg=3)
            cp_aic(j1, j2) = ba
            cp_bic(j1, j2) = bb

            if (print_param) call print_estimated_parameters_dates(max_m, parent, &
                df_ret%values(2:, j1), df_ret%values(2:, j2), date_labels)
        end do
    end do

    ! ── summary matrices ─────────────────────────────────────────────────────────
    ! Diagonal = 0 (not applicable as a pairwise quantity).
    ! Above diagonal = BIC changepoints; below diagonal = AIC changepoints.
    ! The diagonal of cp_aic/cp_bic holds the variance results; we print those
    ! separately since the matrix convention reserves the diagonal for zero.

    print "(/,a)", "Summary: # variance/covariance changepoints chosen by BIC " // &
        "(above diagonal) and AIC (below diagonal)"
    print "(a)", "  Diagonal: variance changepoints (single value if AIC=BIC, else AIC/BIC)"
    print *

    ! Header
    write(*, "(a8)", advance='no') ""
    do j2 = 1, n_col
        write(*, "(a8)", advance='no') trim(df_ret%columns(j2))
    end do
    print *

    do j1 = 1, n_col
        write(*, "(a8)", advance='no') trim(df_ret%columns(j1))
        do j2 = 1, n_col
            if (j1 == j2) then
                if (cp_aic(j1,j1) == cp_bic(j1,j1)) then
                    write(*, "(i8)", advance='no') cp_aic(j1, j1)
                else
                    write(*, "(i4,'/',i2)", advance='no') cp_aic(j1, j1), cp_bic(j1, j1)
                end if
            else if (j2 > j1) then
                write(*, "(i8)", advance='no') cp_bic(j1, j2)   ! above diagonal: BIC
            else
                write(*, "(i8)", advance='no') cp_aic(j2, j1)   ! below diagonal: AIC
            end if
        end do
        print *
    end do

    deallocate(cost, dp_table, parent, z, cp_aic, cp_bic)
    call print_wall_time(t_start)
end program xreturns_cov_cp
