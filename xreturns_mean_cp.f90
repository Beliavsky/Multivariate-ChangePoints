!> program: xreturns_mean_cp
!!
!! reads asset prices, computes returns, and for each asset finds changepoints
!! in the return mean using a profile normal mean-shift model applied to the
!! return series z_t = r_t.
!!
!! cost function per segment: m/2 * log(sample_variance(z))
!! this profiles over both the mean and variance of z within each segment,
!! so it is mainly a mean-shift detector, but it can also react to variance shifts.
!!
!! outputs: for each asset, a table of m, ll, aic, bic, changepoint positions.
!! optionally prints estimated segment means and standard deviations by date.
!! at the end: a summary table of number of changepoints chosen by aic and bic.

program xreturns_mean_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: dataframe_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time, sort_int
    use io_utils_mod, only: print_model_selection
    implicit none

    type(dataframe_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp + 1, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv", &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    logical, parameter :: print_param = .true.
    integer, parameter :: max_days = 260
    logical, parameter :: latest = .true.

    integer :: n, n_col, j, ba, bb
    real(kind=dp), allocatable :: dp_table(:,:), cost(:,:), z(:)
    integer, allocatable :: parent(:,:), cp_aic(:), cp_bic(:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: date_labels(:)

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

    n = size(df_ret%values, 1) - 1
    n_col = ncol(df_ret)
    date_labels = df_ret%index%to_str()

    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n), z(n))
    allocate(cp_aic(n_col), cp_bic(n_col))
    cp_aic = 0
    cp_bic = 0

    do j = 1, n_col
        z = df_ret%values(2:, j)

        print "(/,a)", "changes in mean of " // trim(df_ret%columns(j))

        cost = mean_shift_cost_matrix(z, min_seg_len=min_seg_len)
        call solve_changepoints(max_m, cost, dp_table, parent)
        call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, params_per_seg=3)
        cp_aic(j) = ba
        cp_bic(j) = bb

        if (print_param) call print_estimated_mean_dates(max_m, parent, z, date_labels)
    end do

    print "(/,a)", "summary: number of mean changepoints chosen by aic and bic"
    print "(a10,2a8)", "series", "aic", "bic"
    do j = 1, n_col
        print "(a10,2i8)", trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
    end do

    deallocate(cost, dp_table, parent, z, cp_aic, cp_bic)
    call print_wall_time(t_start)

contains

    !> print estimated segment means and standard deviations using date labels
    subroutine print_estimated_mean_dates(max_m, parent, x, dates)
        integer, intent(in) :: max_m
        integer, intent(in) :: parent(:, :)
        real(kind=dp), intent(in) :: x(:)
        character(len=*), intent(in) :: dates(:)

        integer :: n, m, cp, k, seg_start, seg_end, nseg
        integer :: cps(max_m)
        real(kind=dp) :: mx, sdx

        n = size(parent, 1)

        do m = 1, max_m
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))

            print "(a,i2)", "estimated parameters for m =", m
            print "(a12,a12,a10,a12,a12)", "start", "end", "#obs", "mean", "sd"

            seg_start = 1
            do k = 1, m
                if (k == m) then
                    seg_end = n
                else
                    seg_end = cps(k)
                end if

                nseg = seg_end - seg_start + 1
                mx = sum(x(seg_start:seg_end)) / real(nseg, dp)
                if (nseg > 1) then
                    sdx = sqrt(sum((x(seg_start:seg_end) - mx)**2) / real(nseg - 1, dp))
                else
                    sdx = 0.0_dp
                end if

                print "(a12,a12,i10,2f12.4)", dates(seg_start), dates(seg_end), nseg, mx, sdx
                seg_start = seg_end + 1
            end do
            print *
        end do
    end subroutine print_estimated_mean_dates

end program xreturns_mean_cp
