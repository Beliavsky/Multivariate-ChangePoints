module io_utils_mod
    use kind_mod, only: dp
    use basic_stats_mod, only: mean
    use util_mod, only: sort_int
    implicit none
    public :: print_true_params, print_model_selection, print_estimated_parameters, print_estimated_parameters_dates

contains

    !> Prints the true simulation parameters and empirical correlations for segments.
    subroutine print_true_params(true_cps, corr_true, x, y)
        integer, intent(in) :: true_cps(:)        ! Changepoint locations
        real(kind=dp), intent(in) :: corr_true(:) ! True correlation per segment
        real(kind=dp), intent(in) :: x(:), y(:)   ! Input series (n)
        integer :: n, k, nseg_true, seg_start_true, seg_end_true
        real(kind=dp) :: r_est

        if (size(x) /= size(y)) &
            error stop "print_true_params: size(x) /= size(y)"
        if (size(corr_true) /= size(true_cps) + 1) &
            error stop "print_true_params: size(corr_true) /= size(true_cps) + 1"

        n = size(x)
        nseg_true = size(corr_true)
        print "(a)", "TRUE PARAMETERS"
        print "(a10,a10,a15,a15)", "Start", "End", "Corr_true", "Corr_sim"
        seg_start_true = 1
        do k = 1, nseg_true
            if (k < nseg_true) then
                seg_end_true = true_cps(k)
            else
                seg_end_true = n
            end if
            r_est = sum((x(seg_start_true:seg_end_true)-mean(x(seg_start_true:seg_end_true)))*(y(seg_start_true:seg_end_true)-mean(y(seg_start_true:seg_end_true)))) / &
                    (sqrt(sum((x(seg_start_true:seg_end_true)-mean(x(seg_start_true:seg_end_true)))**2) * sum((y(seg_start_true:seg_end_true)-mean(y(seg_start_true:seg_end_true)))**2)))
            print "(i10,i10,f15.4,f15.4)", seg_start_true, seg_end_true, corr_true(k), r_est
            seg_start_true = seg_end_true + 1
        end do
    end subroutine print_true_params

    !> Prints AIC/BIC statistics for different numbers of segments.
    subroutine print_model_selection(dp_table, parent, best_aic_cp, best_bic_cp, params_per_seg)
        real(kind=dp), intent(in) :: dp_table(:, :) ! Dynamic programming table of costs (n, max_m)
        integer, intent(in) :: parent(:, :)        ! DP parent pointers (n, max_m)
        integer, intent(out), optional :: best_aic_cp  ! # changepoints chosen by AIC
        integer, intent(out), optional :: best_bic_cp  ! # changepoints chosen by BIC
        integer, intent(in), optional :: params_per_seg ! parameters per segment (default 2)
        ! params_per_seg=2: correlation model  → k = 2*m-1 (1 rho/seg + 1 location/break)
        ! params_per_seg=3: mean-shift model   → k = 3*m-1 (mu+sigma^2/seg + 1 location/break)
        integer :: n, max_m, m, cp, k, k_params, pps
        real(kind=dp) :: ll, aic(size(dp_table, 2)), bic(size(dp_table, 2))
        integer :: cps(size(dp_table, 2))
        integer :: best_m_aic, best_m_bic
        real(kind=dp) :: min_aic, min_bic
        pps = 2
        if (present(params_per_seg)) pps = params_per_seg

        if (any(shape(dp_table) /= shape(parent))) &
            error stop "print_model_selection: dp_table and parent shapes are incompatible"

        n = size(dp_table, 1)
        max_m = size(dp_table, 2)

        print "(/,a)", "M      LL         AIC         BIC    ChangePoints"
        do m = 1, max_m
            if (dp_table(n, m) >= 1.0e19_dp) then
                aic(m) = 1.0e20_dp
                bic(m) = 1.0e20_dp
                cycle
            end if
            ll = -dp_table(n, m)
            k_params = pps * m - 1
            aic(m) = -2.0_dp*ll + 2.0_dp * k_params
            bic(m) = -2.0_dp*ll + real(k_params, dp) * log(real(n, dp))

            write(*, "(i2, 3f12.2, 4x)", advance='no') m, ll, aic(m), bic(m)
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))
            do k = 1, m-1
                write(*, "(i0, ' ')", advance='no') cps(k)
            end do
            print *
        end do
        
        min_aic = 1.0e20_dp
        min_bic = 1.0e20_dp
        best_m_aic = 1
        best_m_bic = 1
        do m = 1, max_m
            if (aic(m) < min_aic) then
                min_aic = aic(m)
                best_m_aic = m
            end if
            if (bic(m) < min_bic) then
                min_bic = bic(m)
                best_m_bic = m
            end if
        end do
        print "(/,a,2(1x,i0),/)", "Changepoints chosen by AIC, BIC:", best_m_aic - 1, best_m_bic - 1
        if (present(best_aic_cp)) best_aic_cp = best_m_aic - 1
        if (present(best_bic_cp)) best_bic_cp = best_m_bic - 1
    end subroutine print_model_selection

    !> Prints the estimated parameters (breakpoints and correlations) for each model using indices.
    subroutine print_estimated_parameters(max_m, parent, x, y)
        integer, intent(in) :: max_m              ! Max changepoints
        integer, intent(in) :: parent(:, :)       ! DP parent pointers
        real(kind=dp), intent(in) :: x(:), y(:)   ! Input series
        integer :: n, m, cp, k, seg_start, seg_end
        integer :: cps(max_m)
        real(kind=dp) :: r_est

        if (size(x) /= size(y)) &
            error stop "print_estimated_parameters: size(x) /= size(y)"

        n = size(parent, 1)
        do m = 1, max_m
            ! Reconstruct and sort
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))

            print "(a,i2)", "estimated parameters for m =", m
            print "(a10,a10,a10)", "Start", "End", "Corr"
            seg_start = 1
            do k = 1, m
                if (k == m) then
                    seg_end = n
                else
                    seg_end = cps(k)
                end if
                r_est = sum((x(seg_start:seg_end)-mean(x(seg_start:seg_end)))*(y(seg_start:seg_end)-mean(y(seg_start:seg_end)))) / &
                        (sqrt(sum((x(seg_start:seg_end)-mean(x(seg_start:seg_end)))**2) * sum((y(seg_start:seg_end)-mean(y(seg_start:seg_end)))**2)))
                print "(i10,i10,f10.4)", seg_start, seg_end, r_est
                seg_start = seg_end + 1
            end do
            print *
        end do
    end subroutine print_estimated_parameters

    !> Prints the estimated parameters (breakpoints and correlations) for each model using date labels.
    subroutine print_estimated_parameters_dates(max_m, parent, x, y, dates)
        integer, intent(in) :: max_m              ! Max changepoints
        integer, intent(in) :: parent(:, :)       ! DP parent pointers
        real(kind=dp), intent(in) :: x(:), y(:)   ! Input series
        character(len=*), intent(in) :: dates(:)  ! Date labels
        integer :: n, m, cp, k, seg_start, seg_end
        integer :: cps(max_m)
        real(kind=dp) :: r_est, mx, my, sdx, sdy, covar
        real(kind=dp), allocatable :: x_seg(:), y_seg(:)

        if (size(x) /= size(y)) &
            error stop "print_estimated_parameters_dates: size(x) /= size(y)"

        n = size(parent, 1)
        allocate(x_seg(n), y_seg(n))
        do m = 1, max_m
            ! Reconstruct and sort
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))

            print "(a,i2)", "estimated parameters for m =", m
            print "(a12,a12,a10,a10,a10,a10,a10,a10,a10)", "Start", "End", "#obs", "Corr", "Covar", "sd_x", "sd_y", "mean_x", "mean_y"
            seg_start = 1
            do k = 1, m
                if (k == m) then
                    seg_end = n
                else
                    seg_end = cps(k)
                end if
                
                ! Extract segment
                x_seg(1:seg_end-seg_start+1) = x(seg_start:seg_end)
                y_seg(1:seg_end-seg_start+1) = y(seg_start:seg_end)
                mx = mean(x_seg(1:seg_end-seg_start+1))
                my = mean(y_seg(1:seg_end-seg_start+1))
                sdx = sqrt(sum((x_seg(1:seg_end-seg_start+1)-mx)**2) / max(1, (seg_end-seg_start)))
                sdy = sqrt(sum((y_seg(1:seg_end-seg_start+1)-my)**2) / max(1, (seg_end-seg_start)))
                covar = sum((x_seg(1:seg_end-seg_start+1)-mx)*(y_seg(1:seg_end-seg_start+1)-my)) / max(1, (seg_end-seg_start))
                r_est = covar / (sdx * sdy)
                
                ! r_est, covar, sdx, sdy, mx, my
                print "(a12,a12,i10,f10.4,5f10.4)", dates(seg_start), dates(seg_end), seg_end - seg_start + 1, r_est, covar, sdx, sdy, mx, my

                seg_start = seg_end + 1
            end do
            print *
        end do
        deallocate(x_seg, y_seg)
    end subroutine print_estimated_parameters_dates
end module io_utils_mod
