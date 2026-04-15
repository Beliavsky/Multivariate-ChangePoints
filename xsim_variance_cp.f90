!> Program: xsim_variance_cp
!!
!! Evaluates two variance-changepoint criteria via Monte Carlo simulation.
!! Generates Gaussian returns with known piecewise-constant variance and
!! compares how well BIC recovers the true changepoints under each criterion:
!!
!!   criterion 1:  z_t = r_t^2        (r-squared)
!!   criterion 2:  z_t = log(c + r_t^2)  (log-variance with offset c)
!!
!! For each criterion reports:
!!   mean # changepoints found (vs. true)
!!   fraction of simulations finding exactly the right number
!!   mean absolute location error when count is exact (obs units)
!!   mean distance from each found CP to its nearest true CP (all sims)
!!   histogram of the count of found changepoints

program xsim_variance_cp
    use kind_mod,        only: dp, long_int
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints, segment_ends
    use util_mod,        only: print_wall_time
    implicit none

    ! ── simulation parameters ────────────────────────────────────────────────────
    integer,  parameter :: n_sim      = 1000    ! Monte Carlo replications
    integer,  parameter :: n_per_seg  = 500     ! obs per segment (equal lengths)
    real(dp), parameter :: sigmas(*)  = [1.0_dp, 1.5_dp, 1.0_dp]  ! sigma per segment
    integer,  parameter :: min_seg_len = 50
    integer,  parameter :: max_cp     = 20
    real(dp), parameter :: var_offset = 0.01_dp  ! c in log(c + r^2)
    logical,  parameter :: write_data = .true.  ! write last simulation's data to file
    character(len=*), parameter :: data_file = "sim_variance_cp.txt"

    ! ── derived constants ────────────────────────────────────────────────────────
    integer, parameter :: n_seg     = size(sigmas)
    integer, parameter :: n_cp_true = n_seg - 1
    integer, parameter :: n         = n_seg * n_per_seg
    integer, parameter :: max_m     = max_cp + 1
    integer, parameter :: hist_hi   = n_cp_true + 4  ! histogram rows 0..hist_hi (last = ">=")

    ! ── arrays ───────────────────────────────────────────────────────────────────
    integer  :: true_cps(n_cp_true)
    real(dp), allocatable :: r(:), z(:), dp_tab(:,:)
    integer,  allocatable :: parent(:,:)
    integer               :: seg_buf(max_m)          ! segment_ends workspace

    ! ── accumulators [1 = r^2,  2 = log(c+r^2)] ─────────────────────────────────
    integer  :: n_found
    integer  :: sum_found(2), cnt_exact(2), cnt_any(2)
    real(dp) :: sum_err_exact(2), sum_dist_any(2)
    integer  :: hist(0:hist_hi, 2)

    integer(kind=long_int) :: t_start
    integer :: isim, iseg, k, ic, i, ki
    real(dp) :: d, dmin

    call system_clock(t_start)

    ! true changepoints: end of each non-final segment
    do iseg = 1, n_cp_true
        true_cps(iseg) = iseg * n_per_seg
    end do

    allocate(r(n), z(n), dp_tab(n, max_m), parent(n, max_m))

    sum_found     = 0;  cnt_exact  = 0;  cnt_any   = 0
    sum_err_exact = 0.0_dp;  sum_dist_any = 0.0_dp
    hist          = 0

    ! ── header ───────────────────────────────────────────────────────────────────
    print "('n_sim=',i0,', n_per_seg=',i0,', n=',i0,', n_cp_true=',i0)", &
        n_sim, n_per_seg, n, n_cp_true
    print "('sigmas:', *(1x, f0.3))", sigmas
    print "('true changepoints:', *(1x, i0))", true_cps
    print "('min_seg_len=',i0,', max_cp=',i0,', var_offset=',g0)", &
        min_seg_len, max_cp, var_offset

    ! ── main simulation loop ─────────────────────────────────────────────────────
    do isim = 1, n_sim

        ! generate returns: r(t) ~ N(0, sigma(seg)^2)
        do iseg = 1, n_seg
            do i = (iseg-1)*n_per_seg + 1, iseg*n_per_seg
                r(i) = sigmas(iseg) * randn()
            end do
        end do

        ! evaluate both criteria
        do ic = 1, 2
            if (ic == 1) then
                z = r**2
            else
                z = log(var_offset + r**2)
            end if

            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(z, min_seg_len=min_seg_len), dp_tab, parent)
            n_found = bic_best(dp_tab, n, pps=3)

            sum_found(ic) = sum_found(ic) + n_found
            hist(min(n_found, hist_hi), ic) = hist(min(n_found, hist_hi), ic) + 1

            if (n_found > 0) then
                cnt_any(ic) = cnt_any(ic) + 1
                seg_buf(1:n_found+1) = segment_ends(parent, n_found+1)
                ! seg_buf(1:n_found) = changepoint positions; seg_buf(n_found+1) = n

                ! mean distance from each found CP to its nearest true CP
                d = 0.0_dp
                do k = 1, n_found
                    dmin = huge(1.0_dp)
                    do ki = 1, n_cp_true
                        dmin = min(dmin, abs(real(seg_buf(k) - true_cps(ki), dp)))
                    end do
                    d = d + dmin
                end do
                sum_dist_any(ic) = sum_dist_any(ic) + d / n_found

                ! exact count: mean |found_k - true_k| matched in sorted order
                if (n_found == n_cp_true) then
                    cnt_exact(ic) = cnt_exact(ic) + 1
                    d = 0.0_dp
                    do k = 1, n_cp_true
                        d = d + abs(real(seg_buf(k) - true_cps(k), dp))
                    end do
                    sum_err_exact(ic) = sum_err_exact(ic) + d / n_cp_true
                end if
            end if
        end do
    end do

    ! ── results ──────────────────────────────────────────────────────────────────
    print "(/,a45, 2a16)", " ", "r^2", "log(c+r^2)"
    print "(a45, 2f16.3)", "Mean # changepoints found:", &
        real(sum_found(1), dp)/n_sim, real(sum_found(2), dp)/n_sim
    print "(a45, 2f16.3)", "Fraction with exact count:", &
        real(cnt_exact(1), dp)/n_sim, real(cnt_exact(2), dp)/n_sim
    print "(a45, 2f16.1)", "Mean location error when exact (obs):", &
        safe_mean(sum_err_exact(1), cnt_exact(1)), &
        safe_mean(sum_err_exact(2), cnt_exact(2))
    print "(a45, 2f16.1)", "Mean dist found->nearest true CP (obs):", &
        safe_mean(sum_dist_any(1), cnt_any(1)), &
        safe_mean(sum_dist_any(2), cnt_any(2))

    print "(/,'Distribution of # changepoints found:')"
    print "(a6, 2a16)", "n_cp", "r^2", "log(c+r^2)"
    do i = 0, hist_hi
        if (i < hist_hi) then
            print "(i6, 2f16.3)", i, &
                real(hist(i,1), dp)/n_sim, real(hist(i,2), dp)/n_sim
        else
            print "(i5,'+ ', 2f16.3)", i, &
                real(hist(i,1), dp)/n_sim, real(hist(i,2), dp)/n_sim
        end if
    end do

    if (write_data) then
        open(newunit=k, file=data_file, status='replace')
        write(k, "(a5, a4, 3a14)") "t", "seg", "r", "r^2", "log(c+r^2)"
        do i = 1, n
            iseg = (i-1)/n_per_seg + 1
            write(k, "(i5, i4, 3f14.6)") i, iseg, r(i), r(i)**2, log(var_offset + r(i)**2)
        end do
        close(k)
        print "('wrote ',i0,' rows to ',a)", n, data_file
    end if

    deallocate(r, z, dp_tab, parent)
    call print_wall_time(t_start)

contains

    function randn() result(x)
        !! Box-Muller normal random variate N(0,1).
        real(dp) :: x, u1, u2
        real(dp), parameter :: twopi = 6.2831853071795864769_dp
        do
            call random_number(u1)
            if (u1 > 0.0_dp) exit
        end do
        call random_number(u2)
        x = sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
    end function randn

    pure function bic_best(dp_tab, nt, pps) result(best_cp)
        !! BIC-optimal number of changepoints from dp_tab (no printing).
        real(dp), intent(in) :: dp_tab(:,:)
        integer,  intent(in) :: nt, pps
        integer :: best_cp, m
        real(dp) :: bic, min_bic
        min_bic = huge(1.0_dp)
        best_cp = 0
        do m = 1, size(dp_tab, 2)
            if (dp_tab(nt, m) >= 1.0e19_dp) cycle
            bic = 2.0_dp * dp_tab(nt, m) + real(pps*m - 1, dp) * log(real(nt, dp))
            if (bic < min_bic) then
                min_bic = bic
                best_cp = m - 1
            end if
        end do
    end function bic_best

    pure function safe_mean(s, cnt) result(x)
        !! s / cnt, or 0 if cnt = 0.
        real(dp), intent(in) :: s
        integer,  intent(in) :: cnt
        real(dp) :: x
        x = merge(s / cnt, 0.0_dp, cnt > 0)
    end function safe_mean

end program xsim_variance_cp
