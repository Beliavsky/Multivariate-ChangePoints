!> Program: xcovar_sim
!!
!! Simulates a multivariate normal time series with piecewise-constant
!! covariance matrix, then uses dynamic programming to find changepoints
!! in the joint covariance matrix (same algorithm as xreturns_covmat_cp).
!!
!! True parameters are printed first; estimated parameters for the
!! BIC-chosen model are printed after model selection.

program xcovar_sim
    use kind_mod, only: dp, long_int
    use changepoint_mod, only: multivar_cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time, sort_int
    use sim_changepoint_mod, only: generate_multivar_series
    use io_utils_mod, only: print_model_selection
    implicit none

    ! ── simulation parameters ─────────────────────────────────────────────────
    integer,  parameter :: p         = 3
    integer,  parameter :: n         = 2000
    integer,  parameter :: nseg_true = 3
    integer,  parameter :: true_cps(nseg_true-1) = [500, 1200]
    integer,  parameter :: max_cp    = 20, max_m = max_cp + 1, min_seg_len = 50

    ! Off-diagonal correlations, lower triangle row-major: (2,1),(3,1),(3,2)
    ! per segment (column).
    real(kind=dp), parameter :: rho_true(p*(p-1)/2, nseg_true) = reshape([ &
         0.30_dp,  0.10_dp,  0.50_dp, &   ! segment 1
         0.70_dp,  0.60_dp,  0.80_dp, &   ! segment 2 — high positive
        -0.40_dp,  0.10_dp, -0.20_dp  ], &! segment 3 — mixed
        shape=[p*(p-1)/2, nseg_true])

    ! Daily standard deviations (fractional, e.g. 0.010 = 1% → ~16% annual)
    real(kind=dp), parameter :: sd_true(p, nseg_true) = reshape([ &
        0.010_dp, 0.012_dp, 0.015_dp, &   ! segment 1
        0.020_dp, 0.018_dp, 0.014_dp, &   ! segment 2 — higher vol
        0.008_dp, 0.016_dp, 0.013_dp  ], &! segment 3 — lower / mixed
        shape=[p, nseg_true])

    character(len=2), parameter :: col_names(p) = ['A ', 'B ', 'C ']
    real(kind=dp),    parameter :: scale_ret = 1.0_dp  ! unscaled simulation
    ! ─────────────────────────────────────────────────────────────────────────

    real(kind=dp) :: cov_true(p, p, nseg_true)
    real(kind=dp), allocatable :: R(:,:), dp_table(:,:), cost(:,:)
    integer, allocatable :: parent(:,:), seg_ends(:)
    integer(kind=long_int) :: t_start
    integer :: i, k, a, b, seg, idx, m_segs, seg_start, best_aic, best_bic, pps, cp_idx
    character(len=10) :: ret_dates(n)

    call system_clock(t_start)

    ! ── build true covariance matrices ────────────────────────────────────────
    do seg = 1, nseg_true
        do a = 1, p
            cov_true(a, a, seg) = sd_true(a, seg)**2
        end do
        do a = 2, p
            do b = 1, a-1
                idx = (a-1)*(a-2)/2 + b
                cov_true(a, b, seg) = rho_true(idx, seg) * sd_true(a, seg) * sd_true(b, seg)
                cov_true(b, a, seg) = cov_true(a, b, seg)
            end do
        end do
    end do

    ! ── synthetic date labels (day index as string) ───────────────────────────
    do i = 1, n
        write(ret_dates(i), "(i0)") i
    end do

    ! ── generate data ─────────────────────────────────────────────────────────
    allocate(R(n, p))
    call generate_multivar_series(true_cps, cov_true, R)

    ! ── print true parameters ─────────────────────────────────────────────────
    print "(a)", "TRUE PARAMETERS"
    seg_start = 1
    do seg = 1, nseg_true
        i = merge(true_cps(seg), n, seg < nseg_true)  ! last obs of segment
        call print_segment(seg, seg_start, i, R, col_names, ret_dates, p, &
            cov_true(:,:,seg))
        seg_start = i + 1
    end do

    ! ── cost matrix, DP, model selection ─────────────────────────────────────
    print "(/,a)", "ESTIMATED PARAMETERS"
    print "('building ',i0,'x',i0,' joint cost matrix...')", n, n
    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n))
    cost = multivar_cost_matrix(R, min_seg_len=min_seg_len)
    call solve_changepoints(max_m, cost, dp_table, parent)
    deallocate(cost)

    pps = p * (p + 3) / 2 + 1
    call print_model_selection(dp_table, parent, &
        best_aic_cp=best_aic, best_bic_cp=best_bic, params_per_seg=pps)

    ! ── estimated segment structure for BIC-chosen model ─────────────────────
    m_segs = best_bic + 1
    allocate(seg_ends(m_segs))
    cp_idx = n
    seg_ends(m_segs) = n
    do k = m_segs, 2, -1
        seg_ends(k-1) = parent(cp_idx, k)
        cp_idx = seg_ends(k-1)
    end do
    if (m_segs > 1) call sort_int(seg_ends(1:m_segs-1))

    print "(/,'Covariance structure for BIC-chosen model (',i0,' segment',a,')')", &
        m_segs, merge('s',' ', m_segs /= 1)
    seg_start = 1
    do k = 1, m_segs
        call print_segment(k, seg_start, seg_ends(k), R, col_names, ret_dates, p)
        seg_start = seg_ends(k) + 1
    end do

    deallocate(dp_table, parent, R, seg_ends)
    call print_wall_time(t_start)

contains

    !> Print correlation matrix and annualised std devs for one segment.
    !! If cov_known is present, print it as the "true" matrix alongside the
    !! empirical one; otherwise print the empirical matrix only.
    subroutine print_segment(k, i0, i1, R, col_names, dates, p, cov_known)
        integer,          intent(in)           :: k, i0, i1, p
        real(kind=dp),    intent(in)           :: R(:,:)
        character(len=*), intent(in)           :: col_names(:), dates(:)
        real(kind=dp),    intent(in), optional :: cov_known(p, p)
        integer       :: a, b, m
        real(kind=dp) :: S(p,p), mu(p), sd(p), r_ab
        real(kind=dp), parameter :: ann = 15.87401_dp   ! sqrt(252)

        m = i1 - i0 + 1
        print "(/,'Segment ',i0,': obs ',a,' to ',a,' (',i0,' obs)')", &
            k, trim(dates(i0)), trim(dates(i1)), m

        ! ── true matrix (if supplied) ─────────────────────────────────────────
        if (present(cov_known)) then
            print "(4x,a)", "True correlation / Ann.StdDev:"
            write(*, "(8x)", advance='no')
            do a = 1, p
                write(*, "(a8)", advance='no') trim(col_names(a))
            end do
            print *
            do a = 1, p
                write(*, "(4x,a4)", advance='no') trim(col_names(a))
                do b = 1, a
                    if (cov_known(a,a) > 0.0_dp .and. cov_known(b,b) > 0.0_dp) then
                        r_ab = cov_known(a,b) / sqrt(cov_known(a,a) * cov_known(b,b))
                    else
                        r_ab = 0.0_dp
                    end if
                    write(*, "(f8.3)", advance='no') r_ab
                end do
                print *
            end do
            write(*, "(4x,a4)", advance='no') '*SD*'
            do a = 1, p
                write(*, "(f8.3)", advance='no') sqrt(cov_known(a,a)) * ann / scale_ret
            end do
            print *
        end if

        ! ── empirical matrix ──────────────────────────────────────────────────
        do a = 1, p
            mu(a) = sum(R(i0:i1, a)) / m
        end do
        do a = 1, p
            do b = a, p
                S(a,b) = sum((R(i0:i1,a)-mu(a)) * (R(i0:i1,b)-mu(b))) / m
                S(b,a) = S(a,b)
            end do
        end do
        do a = 1, p
            sd(a) = sqrt(max(S(a,a), 0.0_dp))
        end do

        print "(4x,a)", "Empirical correlation / Ann.StdDev:"
        write(*, "(8x)", advance='no')
        do a = 1, p
            write(*, "(a8)", advance='no') trim(col_names(a))
        end do
        print *
        do a = 1, p
            write(*, "(4x,a4)", advance='no') trim(col_names(a))
            do b = 1, a
                if (sd(a) > 0.0_dp .and. sd(b) > 0.0_dp) then
                    r_ab = S(a,b) / (sd(a) * sd(b))
                else
                    r_ab = 0.0_dp
                end if
                write(*, "(f8.3)", advance='no') r_ab
            end do
            print *
        end do
        write(*, "(4x,a4)", advance='no') '*SD*'
        do a = 1, p
            write(*, "(f8.3)", advance='no') sd(a) * ann / scale_ret
        end do
        print *
    end subroutine print_segment

end program xcovar_sim
