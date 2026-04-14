!> program: xreturns_variance_pwl_fast
!!
!! reads asset prices, computes returns, and for each asset fits a continuous
!! piecewise linear variance curve to squared returns.
!!
!! model:
!!   z_t = r_t^2
!!   v_t is continuous piecewise linear in t, with knots at changepoints
!!   z_t is fit by least squares to v_t using a linear spline basis
!!
!! unlike xreturns_variance_cp.f90, this is not an exact dynamic-programming
!! changepoint solver. continuity of the variance curve couples adjacent
!! segments, so this program uses greedy knot insertion plus local refinement.
!!
!! this is the faster heuristic companion to xreturns_variance_pwl.f90.
!!
!! model selection is based on the gaussian spline-fit loglikelihood
!! for squared returns around the fitted variance curve.
!!
!! outputs: for each asset, a table of m, ll, aic, bic, changepoint positions.
!! optionally prints fitted knot variances by date.

program xreturns_variance_pwl_fast
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: dataframe_index_date, nrow, ncol, operator(*)
    use util_mod, only: print_wall_time, sort_int
    implicit none

    type(dataframe_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, min_seg_len = 50, max_assets = 1000
    integer, parameter :: max_refine_iter = 3
    character(len=*), parameter :: prices_file = 'spy_efa_eem_tlt.csv', &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    real(kind=dp), parameter :: min_var = 1.0e-8_dp
    logical, parameter :: print_param = .false.
    integer, parameter :: max_days = 260
    logical, parameter :: latest = .true.

    integer :: n, n_col, j, max_cp_eff, best_aic_m, best_bic_m, n_models_fitted
    real(kind=dp), allocatable :: r(:), ll(:), aic(:), bic(:)
    integer, allocatable :: cp_store(:, :), cp_aic(:), cp_bic(:)
    real(kind=dp), allocatable :: fitted_best(:, :)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: date_labels(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print '(*(1x,a,1x,i0))', 'read', nrow(df_px), 'days and', ncol(df_px), &
        'columns from ' // trim(prices_file)

    if (max_days > 0) then
        df_px = df_px%keep_rows(max_days + 1, latest=latest)
        print '(''using '',a,'' '',i0,'' returns ('',a,'' to '',a,'')'')', &
            merge('latest  ','earliest', latest), nrow(df_px) - 1, &
            trim(df_px%index(2)%to_str()), trim(df_px%index(nrow(df_px))%to_str())
    end if

    df_ret = scale_ret * df_px%pct_change()
    if (scale_ret /= 1.0_dp) print '(''return scaling: '', f0.4)', scale_ret

    n = size(df_ret%values, 1) - 1
    n_col = ncol(df_ret)
    max_cp_eff = min(max_cp, max(0, n / min_seg_len - 1))

    if (max_cp_eff < max_cp) then
        print '(/,a,i0,a,i0)', 'effective max_cp reduced from ', max_cp, ' to ', max_cp_eff
    end if

    allocate(r(n), ll(max_cp_eff + 1), aic(max_cp_eff + 1), bic(max_cp_eff + 1))
    allocate(cp_store(max_cp_eff, max_cp_eff + 1), cp_aic(n_col), cp_bic(n_col))
    allocate(fitted_best(n, max_cp_eff + 1))
    cp_store = 0
    cp_aic = 0
    cp_bic = 0
    fitted_best = 0.0_dp

    if (print_param) then
        allocate(date_labels(n))
        do j = 1, n
            date_labels(j) = df_px%index(j + 1)%to_str()
        end do
    end if

    do j = 1, n_col
        r = df_ret%values(2:, j)

        print '(/,a)', 'continuous piecewise linear variance of ' // trim(df_ret%columns(j))

        call fit_series_pwl(r, min_seg_len, max_cp_eff, max_refine_iter, &
            ll, aic, bic, cp_store, fitted_best, n_models_fitted)

        call print_model_selection_pwl(ll, aic, bic, cp_store, n_models_fitted, best_aic_m, best_bic_m)
        cp_aic(j) = best_aic_m - 1
        cp_bic(j) = best_bic_m - 1

        if (print_param) then
            call print_knot_values(best_bic_m - 1, cp_store(:, best_bic_m), fitted_best(:, best_bic_m), date_labels)
        end if
    end do

    print '(/,a)', 'summary: # variance changepoints chosen by aic and bic'
    print '(a10,2a8)', 'series', 'aic', 'bic'
    do j = 1, n_col
        print '(a10,2i8)', trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
    end do

    if (allocated(date_labels)) deallocate(date_labels)
    deallocate(r, ll, aic, bic, cp_store, cp_aic, cp_bic, fitted_best)
    call print_wall_time(t_start)

contains

    subroutine fit_series_pwl(r, min_seg_len, max_cp_eff, max_refine_iter, ll, aic, bic, cp_store, fitted_store, n_models_fitted)
    ! fit a sequence of continuous piecewise linear variance models for one series.
        real(kind=dp), intent(in) :: r(:)
        integer, intent(in) :: min_seg_len, max_cp_eff, max_refine_iter
        real(kind=dp), intent(out) :: ll(:), aic(:), bic(:)
        integer, intent(out) :: cp_store(:, :)
        real(kind=dp), intent(out) :: fitted_store(:, :)
        integer, intent(out) :: n_models_fitted

        integer :: n, m, ncp
        real(kind=dp), allocatable :: y(:), fitted(:)
        integer, allocatable :: cps(:)

        n = size(r)
        allocate(y(n), fitted(n), cps(max_cp_eff))
        y = r**2
        cps = 0
        cp_store = 0
        fitted_store = 0.0_dp

        ncp = 0
        call fit_y_to_cps(y, cps(1:ncp), fitted)
        call score_model(y, fitted, ncp, ll(1), aic(1), bic(1))
        fitted_store(:, 1) = fitted

        n_models_fitted = 1
        do m = 2, max_cp_eff + 1
            ncp = m - 1
            call add_best_cp(y, cps(1:ncp-1), min_seg_len, cps(ncp))
            if (cps(ncp) < 0) exit
            call sort_int(cps(1:ncp))
            call refine_cps(y, cps(1:ncp), min_seg_len, max_refine_iter)
            call fit_y_to_cps(y, cps(1:ncp), fitted)
            call score_model(y, fitted, ncp, ll(m), aic(m), bic(m))
            cp_store(1:ncp, m) = cps(1:ncp)
            fitted_store(:, m) = fitted
            n_models_fitted = m
        end do

        if (n_models_fitted < max_cp_eff + 1) then
            ll(n_models_fitted+1:max_cp_eff+1) = -huge(1.0_dp)
            aic(n_models_fitted+1:max_cp_eff+1) = huge(1.0_dp)
            bic(n_models_fitted+1:max_cp_eff+1) = huge(1.0_dp)
        end if

        deallocate(y, fitted, cps)
    end subroutine fit_series_pwl

    subroutine score_model(y, fitted, ncp, ll, aic, bic)
    ! compute gaussian spline-fit loglikelihood, aic, and bic for squared returns.
        real(kind=dp), intent(in) :: y(:), fitted(:)
        integer, intent(in) :: ncp
        real(kind=dp), intent(out) :: ll, aic, bic
        integer :: n, k_params
        real(kind=dp) :: rss, s2
        real(kind=dp), parameter :: twopi = 6.2831853071795864769_dp

        n = size(y)
        rss = sum((y - fitted)**2)
        s2 = max(rss / real(n, dp), 1.0e-12_dp)
        ll = -0.5_dp * real(n, dp) * (log(twopi * s2) + 1.0_dp)
        k_params = 2 * ncp + 3
        aic = -2.0_dp * ll + 2.0_dp * real(k_params, dp)
        bic = -2.0_dp * ll + log(real(n, dp)) * real(k_params, dp)
    end subroutine score_model

    subroutine add_best_cp(y, cps_old, min_seg_len, best_cp)
    ! add the cp that gives the best least-squares continuous spline fit.
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: cps_old(:), min_seg_len
        integer, intent(out) :: best_cp

        integer :: n, ncp_old, cand
        integer, allocatable :: cps_try(:)
        real(kind=dp) :: rss_best, rss_try
        real(kind=dp), allocatable :: fitted(:)

        n = size(y)
        ncp_old = size(cps_old)
        allocate(cps_try(ncp_old + 1), fitted(n))

        rss_best = huge(1.0_dp)
        best_cp = -1
        do cand = min_seg_len, n - min_seg_len
            cps_try(1:ncp_old) = cps_old
            cps_try(ncp_old + 1) = cand
            call sort_int(cps_try)
            if (.not. valid_cps(cps_try, n, min_seg_len)) cycle
            call fit_y_to_cps(y, cps_try, fitted, rss_try)
            if (rss_try < rss_best) then
                rss_best = rss_try
                best_cp = cand
            end if
        end do

        deallocate(cps_try, fitted)
    end subroutine add_best_cp

    subroutine refine_cps(y, cps, min_seg_len, max_refine_iter)
    ! locally refine changepoints by coordinate search on rss.
        real(kind=dp), intent(in) :: y(:)
        integer, intent(inout) :: cps(:)
        integer, intent(in) :: min_seg_len, max_refine_iter

        integer :: n, ncp, iter, k, cand, lo, hi, best_cand
        real(kind=dp) :: rss_best, rss_try
        real(kind=dp), allocatable :: fitted(:)
        integer, allocatable :: cps_try(:)
        logical :: changed

        n = size(y)
        ncp = size(cps)
        if (ncp == 0) return

        allocate(fitted(n), cps_try(ncp))

        do iter = 1, max_refine_iter
            changed = .false.
            do k = 1, ncp
                cps_try = cps
                call fit_y_to_cps(y, cps_try, fitted, rss_best)
                best_cand = cps(k)

                if (k == 1) then
                    lo = min_seg_len
                else
                    lo = cps(k-1) + min_seg_len
                end if
                if (k == ncp) then
                    hi = n - min_seg_len
                else
                    hi = cps(k+1) - min_seg_len
                end if

                do cand = lo, hi
                    if (cand == cps(k)) cycle
                    cps_try = cps
                    cps_try(k) = cand
                    call fit_y_to_cps(y, cps_try, fitted, rss_try)
                    if (rss_try < rss_best) then
                        rss_best = rss_try
                        best_cand = cand
                    end if
                end do

                if (best_cand /= cps(k)) then
                    cps(k) = best_cand
                    changed = .true.
                end if
            end do
            if (.not. changed) exit
        end do

        call sort_int(cps)
        if (.not. valid_cps(cps, n, min_seg_len)) error stop 'refine_cps: invalid changepoint set'

        deallocate(fitted, cps_try)
    end subroutine refine_cps

    logical function valid_cps(cps, n, min_seg_len)
    ! test whether cps satisfy the minimum segment-length constraint.
        integer, intent(in) :: cps(:), n, min_seg_len
        integer :: k

        valid_cps = .true.
        if (size(cps) == 0) return

        if (cps(1) < min_seg_len) then
            valid_cps = .false.
            return
        end if
        if (n - cps(size(cps)) < min_seg_len) then
            valid_cps = .false.
            return
        end if
        do k = 2, size(cps)
            if (cps(k) - cps(k-1) < min_seg_len) then
                valid_cps = .false.
                return
            end if
        end do
    end function valid_cps

    subroutine fit_y_to_cps(y, cps, fitted, rss)
    ! fit a continuous piecewise linear spline to y for a fixed cp set.
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(out) :: fitted(:)
        real(kind=dp), intent(out), optional :: rss

        integer :: n, p, i, k
        real(kind=dp), allocatable :: xtx(:, :), xty(:), beta(:)
        real(kind=dp) :: t

        n = size(y)
        p = size(cps) + 2
        allocate(xtx(p, p), xty(p), beta(p))
        xtx = 0.0_dp
        xty = 0.0_dp

        do i = 1, n
            t = real(i, dp)
            call accumulate_normal_equations(t, y(i), cps, xtx, xty)
        end do

        call solve_linear_system(xtx, xty, beta)

        do i = 1, n
            t = real(i, dp)
            fitted(i) = beta(1) + beta(2) * t
            do k = 1, size(cps)
                fitted(i) = fitted(i) + beta(k + 2) * max(0.0_dp, t - real(cps(k), dp))
            end do
        end do

        if (present(rss)) rss = sum((y - fitted)**2)

        deallocate(xtx, xty, beta)
    end subroutine fit_y_to_cps

    subroutine accumulate_normal_equations(t, yval, cps, xtx, xty)
    ! add one observation to the normal equations for the spline basis.
        real(kind=dp), intent(in) :: t, yval
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(inout) :: xtx(:, :), xty(:)

        integer :: p, i, j
        real(kind=dp), allocatable :: x(:)

        p = size(cps) + 2
        allocate(x(p))
        x(1) = 1.0_dp
        x(2) = t
        do i = 1, size(cps)
            x(i + 2) = max(0.0_dp, t - real(cps(i), dp))
        end do

        do i = 1, p
            xty(i) = xty(i) + x(i) * yval
            do j = 1, p
                xtx(i, j) = xtx(i, j) + x(i) * x(j)
            end do
        end do

        deallocate(x)
    end subroutine accumulate_normal_equations

    subroutine solve_linear_system(a, b, x)
    ! solve a*x=b by gaussian elimination with partial pivoting.
        real(kind=dp), intent(in) :: a(:, :), b(:)
        real(kind=dp), intent(out) :: x(:)

        integer :: n, i, k, ipiv
        real(kind=dp), allocatable :: aa(:, :), bb(:), rowtmp(:)
        real(kind=dp) :: piv, factor, best

        n = size(b)
        allocate(aa(n, n), bb(n), rowtmp(n))
        aa = a
        bb = b

        do k = 1, n - 1
            ipiv = k
            best = abs(aa(k, k))
            do i = k + 1, n
                if (abs(aa(i, k)) > best) then
                    best = abs(aa(i, k))
                    ipiv = i
                end if
            end do
            if (best <= 1.0e-12_dp) error stop 'solve_linear_system: singular matrix'

            if (ipiv /= k) then
                rowtmp = aa(k, :)
                aa(k, :) = aa(ipiv, :)
                aa(ipiv, :) = rowtmp
                piv = bb(k)
                bb(k) = bb(ipiv)
                bb(ipiv) = piv
            end if

            do i = k + 1, n
                factor = aa(i, k) / aa(k, k)
                aa(i, k:n) = aa(i, k:n) - factor * aa(k, k:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do
        if (abs(aa(n, n)) <= 1.0e-12_dp) error stop 'solve_linear_system: singular matrix'

        x(n) = bb(n) / aa(n, n)
        do i = n - 1, 1, -1
            x(i) = (bb(i) - sum(aa(i, i+1:n) * x(i+1:n))) / aa(i, i)
        end do

        deallocate(aa, bb, rowtmp)
    end subroutine solve_linear_system

    subroutine print_model_selection_pwl(ll, aic, bic, cp_store, n_models_fitted, best_aic_m, best_bic_m)
    ! print the model-selection table for the continuous spline fit.
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
    ! print fitted variance values at the start, knots, and end.
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

end program xreturns_variance_pwl_fast
