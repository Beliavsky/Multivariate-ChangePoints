!> program: xreturns_variance_pwl
!!
!! reads asset prices, computes returns, and for each asset fits the exact
!! optimal continuous piecewise linear least-squares model for squared returns.
!!
!! model:
!!   y_t = r_t^2
!!   v_t is continuous piecewise linear in t, with changepoints at internal knots
!!   the fitted path minimizes rss = sum_t (y_t - v_t)^2 exactly for each fixed
!!   number of segments, subject to the minimum segment-length constraint.
!!
!! this is an exact optimizer for the continuous piecewise linear rss objective.
!! it is not the gaussian return loglikelihood for a heteroskedastic variance path.
!!
!! outputs: for each asset, a table of m, ll, aic, bic, changepoint positions.
!! optionally prints fitted variance values at the start, changepoints, and end.

program xreturns_variance_pwl
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: dataframe_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: solve_continuous_pwl_sse, fit_continuous_pwl_given_cps
    use util_mod, only: print_wall_time
    implicit none

    type(dataframe_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = 'spy_efa_eem_tlt.csv', &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    real(kind=dp), parameter :: min_var = 1.0e-8_dp
    logical, parameter :: print_param  = .false.
    integer, parameter :: max_days     = 260
    logical, parameter :: latest       = .true.
    logical, parameter :: resample_ret = .false.  ! .true. = shuffle rows (null hypothesis check)

    integer :: n, n_col, j, max_cp_eff, max_m, best_aic_m, best_bic_m, n_models_fitted
    real(kind=dp), allocatable :: r(:), y(:), sse(:), ll(:), aic(:), bic(:), fitted(:)
    integer, allocatable :: cp_store(:, :), cp_aic(:), cp_bic(:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: date_labels(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print '(*(1x,a,1x,i0))', 'read', nrow(df_px), 'days and', ncol(df_px), &
        'columns from ' // trim(prices_file)

    if (max_days > 0) then
        df_px = df_px%keep_rows(max_days + 1, latest=latest)
    end if

    df_ret = scale_ret * df_px%pct_change(dropna=.true.)
    if (scale_ret /= 1.0_dp) print '(''return scaling: '', f0.4)', scale_ret
    if (resample_ret) then
        df_ret = df_ret%resample()
        print *, "resampled returns!"
    end if

    n = size(df_ret%values, 1)
    n_col = ncol(df_ret)
    max_cp_eff = min(max_cp, max(0, n / min_seg_len - 1))
    max_m = max_cp_eff + 1

    if (max_cp_eff < max_cp) then
        print '(/,a,i0,a,i0)', 'effective max_cp reduced from ', max_cp, ' to ', max_cp_eff
    end if

    allocate(r(n), y(n), fitted(n), sse(max_m), ll(max_m), aic(max_m), bic(max_m))
    allocate(cp_store(max(1, max_m - 1), max_m), cp_aic(n_col), cp_bic(n_col))
    cp_aic = 0
    cp_bic = 0

    if (print_param) then
        allocate(date_labels(n))
        do j = 1, n
            date_labels(j) = df_ret%index(j)%to_str()
        end do
    end if

    do j = 1, n_col
        r = df_ret%values(:, j)
        y = r**2

        print '(/,a)', 'continuous piecewise linear variance of ' // trim(df_ret%columns(j))

        call solve_continuous_pwl_sse(y, max_m, min_seg_len, sse, cp_store, n_models_fitted)
        call score_models_from_sse(y, sse, n_models_fitted, ll, aic, bic)
        call print_model_selection_pwl(ll, aic, bic, cp_store, n_models_fitted, best_aic_m, best_bic_m)
        cp_aic(j) = best_aic_m - 1
        cp_bic(j) = best_bic_m - 1

        if (print_param) then
            call fit_continuous_pwl_given_cps(y, cp_store(1:best_bic_m-1, best_bic_m), fitted)
            call print_knot_values(best_bic_m - 1, cp_store(1:best_bic_m-1, best_bic_m), fitted, date_labels)
        end if
    end do

    print '(/,a)', 'summary: # variance changepoints chosen by aic and bic'
    print '(a10,2a8)', 'series', 'aic', 'bic'
    do j = 1, n_col
        print '(a10,2i8)', trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
    end do

    if (allocated(date_labels)) deallocate(date_labels)
    deallocate(r, y, fitted, sse, ll, aic, bic, cp_store, cp_aic, cp_bic)
    call print_wall_time(t_start)

contains

    subroutine score_models_from_sse(y, sse, n_models_fitted, ll, aic, bic)
    ! compute gaussian spline-fit loglikelihood, aic, and bic from rss values.
        real(kind=dp), intent(in) :: y(:), sse(:)
        integer, intent(in) :: n_models_fitted
        real(kind=dp), intent(out) :: ll(:), aic(:), bic(:)
        integer :: n, m, ncp, k_params
        real(kind=dp) :: sigma2
        real(kind=dp), parameter :: twopi = 6.2831853071795864769_dp

        n = size(y)
        ll = -huge(1.0_dp)
        aic = huge(1.0_dp)
        bic = huge(1.0_dp)

        do m = 1, n_models_fitted
            ncp = m - 1
            sigma2 = max(sse(m) / real(n, dp), 1.0e-12_dp)
            ll(m) = -0.5_dp * real(n, dp) * (log(twopi * sigma2) + 1.0_dp)
            k_params = 2 * ncp + 3
            aic(m) = -2.0_dp * ll(m) + 2.0_dp * real(k_params, dp)
            bic(m) = -2.0_dp * ll(m) + log(real(n, dp)) * real(k_params, dp)
        end do
    end subroutine score_models_from_sse

    subroutine print_model_selection_pwl(ll, aic, bic, cp_store, n_models_fitted, best_aic_m, best_bic_m)
    ! print the model-selection table for the exact continuous spline fit.
        real(kind=dp), intent(in) :: ll(:), aic(:), bic(:)
        integer, intent(in) :: cp_store(:, :), n_models_fitted
        integer, intent(out) :: best_aic_m, best_bic_m

        integer :: m, k
        real(kind=dp) :: min_aic, min_bic

        best_aic_m = 1
        best_bic_m = 1
        min_aic = aic(1)
        min_bic = bic(1)

        print '(/,a)', 'm      ll         aic         bic    changepoints'
        do m = 1, n_models_fitted
            write(*, '(i2,3f12.2,4x)', advance='no') m, ll(m), aic(m), bic(m)
            do k = 1, m - 1
                write(*, '(i0,1x)', advance='no') cp_store(k, m)
            end do
            print *
            if (aic(m) < min_aic) then
                min_aic = aic(m)
                best_aic_m = m
            end if
            if (bic(m) < min_bic) then
                min_bic = bic(m)
                best_bic_m = m
            end if
        end do

        print '(/,a,2(1x,i0),/)', 'changepoints chosen by aic, bic:', best_aic_m - 1, best_bic_m - 1
    end subroutine print_model_selection_pwl

    subroutine print_knot_values(ncp, cps, fitted, dates)
    ! print fitted variance values at the start, changepoints, and end.
        integer, intent(in) :: ncp
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(in) :: fitted(:)
        character(len=*), intent(in) :: dates(:)

        integer :: k, idx

        print '(a)', 'fitted variance at start, changepoints, and end'
        print '(a12,a12)', 'date', 'variance'
        print '(a12,f12.4)', dates(1), max(fitted(1), min_var)
        do k = 1, ncp
            idx = cps(k)
            print '(a12,f12.4)', dates(idx), max(fitted(idx), min_var)
        end do
        print '(a12,f12.4)', dates(size(fitted)), max(fitted(size(fitted)), min_var)
        print *
    end subroutine print_knot_values

end program xreturns_variance_pwl
