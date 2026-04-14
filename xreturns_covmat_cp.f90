!> Program: xreturns_covmat_cp
!!
!! Reads asset prices, computes returns, and finds changepoints in the FULL
!! p×p covariance matrix jointly using dynamic programming.
!!
!! Model: within each segment the returns are i.i.d. multivariate normal with
!! covariance matrix Sigma.  The profile log-likelihood (profiling over both
!! the mean vector and Sigma) gives the cost function
!!
!!     cost(i..j) = m/2 * log|Sigma_hat(i:j)|
!!
!! where Sigma_hat is the biased sample covariance matrix of the segment.
!! This is the multivariate generalisation of the pairwise covariance cost.
!!
!! A single set of changepoints is found for all assets simultaneously,
!! rather than independent per-pair results as in xreturns_cov_cp.f90.
!!
!! The BIC parameter count per segment is p*(p+3)/2 (p means + p*(p+1)/2
!! covariance entries) plus one per changepoint location, giving
!!     k = (p*(p+3)/2 + 1)*m - 1
!! i.e. params_per_seg = p*(p+3)/2 + 1.
!!
!! After model selection the estimated correlation matrix and annualised
!! standard deviations are printed for each segment of the BIC-chosen model.

program xreturns_covmat_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: multivar_cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time, sort_int
    use io_utils_mod, only: print_model_selection
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp + 1, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv", &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    integer, parameter :: max_days = 0     ! 0 = use all returns
    logical, parameter :: latest  = .true. ! if max_days>0: .true.=latest, .false.=earliest

    integer :: n, n_col, k, pps, best_aic, best_bic, m_segs, seg_start, cp_idx
    integer, allocatable :: seg_ends(:)
    real(kind=dp), allocatable :: dp_table(:,:), cost(:,:), R(:,:)
    integer, allocatable :: parent(:,:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: date_labels(:), ret_dates(:)

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

    ! ret_dates(i) is the date of return i (1-indexed), skipping the NaN first row
    date_labels = df_ret%index%to_str()
    allocate(ret_dates(n))
    ret_dates = date_labels(2:n+1)

    ! R(i, a) = return of asset a on day i
    allocate(R(n, n_col))
    R = df_ret%values(2:, :)

    ! ── cost matrix and DP ───────────────────────────────────────────────────────
    print "(/,'building ',i0,'x',i0,' joint cost matrix for ',i0,' assets...')", n, n, n_col
    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n))
    cost = multivar_cost_matrix(R, min_seg_len=min_seg_len)

    print "(a)", "running DP..."
    call solve_changepoints(max_m, cost, dp_table, parent)
    deallocate(cost)

    ! ── model selection ──────────────────────────────────────────────────────────
    ! k = pps*m - 1,  pps = p*(p+3)/2 + 1
    pps = n_col * (n_col + 3) / 2 + 1
    print "('params_per_seg = ',i0,'  (p=',i0,', p*(p+3)/2+1 = ',i0,')')", pps, n_col, pps
    call print_model_selection(dp_table, parent, &
        best_aic_cp=best_aic, best_bic_cp=best_bic, params_per_seg=pps)

    ! ── segment details for BIC-chosen model ─────────────────────────────────────
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
        call print_segment(k, seg_start, seg_ends(k), R, df_ret%columns, ret_dates, n_col)
        seg_start = seg_ends(k) + 1
    end do

    deallocate(dp_table, parent, R, seg_ends)
    call print_wall_time(t_start)

contains

    subroutine print_segment(k, i0, i1, R, col_names, dates, p)
        integer,          intent(in) :: k, i0, i1, p
        real(kind=dp),    intent(in) :: R(:,:)
        character(len=*), intent(in) :: col_names(:), dates(:)
        integer       :: a, b, m
        real(kind=dp) :: S(p,p), mu(p), sd(p), r_ab
        real(kind=dp), parameter :: ann = 15.87401_dp  ! sqrt(252), daily → annual

        m = i1 - i0 + 1
        print "(/,'Segment ',i0,': ',a,' to ',a,' (',i0,' obs)')", &
            k, trim(dates(i0)), trim(dates(i1)), m

        ! Biased sample covariance matrix
        do a = 1, p
            mu(a) = sum(R(i0:i1, a)) / m
        end do
        do a = 1, p
            do b = a, p
                S(a,b) = sum((R(i0:i1,a) - mu(a)) * (R(i0:i1,b) - mu(b))) / m
                S(b,a) = S(a,b)
            end do
        end do
        do a = 1, p
            sd(a) = sqrt(max(S(a,a), 0.0_dp))
        end do

        ! Correlation matrix (lower triangle)
        print "(a)", "  Correlation:"
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

        ! Annualised std devs (unscaled) — aligned with correlation columns above
        write(*, "(4x,a4)", advance='no') '*SD*'
        do a = 1, p
            write(*, "(f8.3)", advance='no') sd(a) * ann / scale_ret
        end do
        print *
    end subroutine print_segment

end program xreturns_covmat_cp
