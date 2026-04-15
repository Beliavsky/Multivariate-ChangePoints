!> Program: xreturns_corrsub_cp
!!
!! Reads asset prices, computes returns, and searches for changepoints using
!! multiple algorithms across all distinct sub-matrices of specified sizes.
!!
!! Toggles:
!!   do_mean     -- changepoints in mean,     per asset (univariate)
!!   do_variance -- changepoints in variance, per asset (univariate)
!!   do_cov      -- changepoints in covariance matrix, for each subset
!!   do_corr     -- changepoints in correlation matrix, for each subset
!!                  (uses EWMA-standardized returns to isolate correlation)
!!
!! sub_sizes specifies which subset sizes to analyze for do_cov / do_corr.
!! All entries must be >= 2.  If C(n_assets, k) > max_subsets only the first
!! max_subsets lexicographic subsets are analyzed and a warning is printed.

program xreturns_corrsub_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: multivar_cost_matrix, mean_shift_cost_matrix, &
        solve_changepoints, segment_ends
    use basic_stats_mod, only: standardize_returns
    use util_mod, only: print_wall_time, next_combination, n_choose_k, join
    use io_utils_mod, only: print_model_selection, keep_obs, corrmat_model_range, &
        print_univar_segments, print_covmat_model, print_corrmat_model
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret
    integer, parameter :: max_cp = 20, max_m = max_cp + 1, min_seg_len = 50, &
        max_assets = 1000
    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv", &
        sym_allowed(*) = [character(len=5) ::]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    integer, parameter  :: max_days = 0      ! 0 = use all data
    logical, parameter  :: latest   = .true. ! if max_days > 0: use latest or earliest
    logical, parameter  :: use_ewma    = .true.
    real(kind=dp), parameter :: ewma_lambda = 0.94_dp
    logical, parameter  :: print_all_models = .false.
    ! print_segs: 0 = BIC-chosen model only
    !             1 = 0 changepoints through BIC-chosen
    !             2 = all models studied
    integer, parameter  :: print_segs = 1
    logical, parameter  :: print_diffs = .true.
    real(kind=dp), parameter :: alpha_diff = 0.05_dp
    integer, parameter  :: max_subsets = 50
    integer, parameter  :: sub_sizes(*) = [2, 5]
    logical, parameter  :: do_mean     = .true.
    logical, parameter  :: do_variance = .true.
    logical, parameter  :: do_cov      = .true.
    logical, parameter  :: do_corr     = .true.

    integer :: n, n_col, j, ks, k, n_total, n_limit, n_sub, pps, m_lo, m_hi, ms, best_bic_sub
    integer, allocatable :: combo(:), cp_aic(:), cp_bic(:), seg_ends(:)
    real(kind=dp), allocatable :: dp_table(:,:), R(:,:), R_std(:,:)
    integer, allocatable :: parent(:,:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: ret_dates(:), date_labels(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read", nrow(df_px), "days and", ncol(df_px), &
        "columns from " // trim(prices_file)
    if (max_days > 0) call keep_obs(df_px, max_days, latest, verbose=.true.)
    df_ret = scale_ret * df_px%pct_change()
    if (scale_ret /= 1.0_dp) print "('return scaling: ', f0.4)", scale_ret

    n     = size(df_ret%values, 1) - 1
    n_col = ncol(df_ret)
    date_labels = df_ret%index%to_str()
    ret_dates   = date_labels(2:n+1)
    R           = df_ret%values(2:, :)

    allocate(dp_table(n, max_m), parent(n, max_m), seg_ends(max_m))

    ! ── print parameters ───────────────────────────────────────────────────────
    print "(/,'do_mean = ',l1,', do_variance = ',l1,', do_cov = ',l1,', do_corr = ',l1)", &
        do_mean, do_variance, do_cov, do_corr
    print "('print_segs = ',i0,'  (0=BIC only, 1=0..BIC, 2=all)')", print_segs
    if (do_cov .or. do_corr) then
        print "('sub_sizes   = ',*(i0,:,' '))", sub_sizes
        print "('max_subsets = ',i0)", max_subsets
    end if
    if (do_corr) print "('use_ewma = ',l1,', ewma_lambda = ',f0.2)", use_ewma, ewma_lambda

    ! ── mean changepoints (one per asset) ──────────────────────────────────────
    if (do_mean) then
        cp_aic = [(0, j=1,n_col)]
        cp_bic = [(0, j=1,n_col)]
        do j = 1, n_col
            print "(/,'changes in mean of ',a)", trim(df_ret%columns(j))
            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(R(:,j), min_seg_len=min_seg_len), &
                dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=cp_aic(j), &
                best_bic_cp=cp_bic(j), params_per_seg=3, print_each=print_all_models)
            call corrmat_model_range(print_segs, cp_bic(j), max_m, m_lo, m_hi)
            do ms = m_lo, m_hi
                if (dp_table(n, ms) >= 1.0e19_dp) cycle
                seg_ends(1:ms) = segment_ends(parent, ms)
                call print_univar_segments(cp_bic(j), seg_ends(1:ms), R(:,j), &
                    trim(df_ret%columns(j)), ret_dates)
            end do
        end do
        print "(/,'summary: # mean changepoints chosen by AIC and BIC')"
        print "(a12,2a8)", "series", "AIC", "BIC"
        do j = 1, n_col
            print "(a12,2i8)", trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
        end do
    end if

    ! ── variance changepoints (one per asset) ──────────────────────────────────
    if (do_variance) then
        cp_aic = [(0, j=1,n_col)]
        cp_bic = [(0, j=1,n_col)]
        do j = 1, n_col
            print "(/,'changes in variance of ',a)", trim(df_ret%columns(j))
            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(R(:,j)**2, min_seg_len=min_seg_len), &
                dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=cp_aic(j), &
                best_bic_cp=cp_bic(j), params_per_seg=3, print_each=print_all_models)
            call corrmat_model_range(print_segs, cp_bic(j), max_m, m_lo, m_hi)
            do ms = m_lo, m_hi
                if (dp_table(n, ms) >= 1.0e19_dp) cycle
                seg_ends(1:ms) = segment_ends(parent, ms)
                call print_univar_segments(cp_bic(j), seg_ends(1:ms), R(:,j)**2, &
                    trim(df_ret%columns(j)), ret_dates)
            end do
        end do
        print "(/,'summary: # variance changepoints chosen by AIC and BIC')"
        print "(a12,2a8)", "series", "AIC", "BIC"
        do j = 1, n_col
            print "(a12,2i8)", trim(df_ret%columns(j)), cp_aic(j), cp_bic(j)
        end do
    end if

    ! ── covariance / correlation changepoints on subsets ───────────────────────
    if (do_cov .or. do_corr) then
        if (do_corr) R_std = standardize_returns(R, use_ewma, ewma_lambda)

        do ks = 1, size(sub_sizes)
            k = sub_sizes(ks)
            if (k < 2) then
                print "('skipping size ',i0,': must be >= 2')", k
                cycle
            end if
            if (k > n_col) then
                print "('skipping size ',i0,': exceeds n_col = ',i0)", k, n_col
                cycle
            end if
            pps     = k * (k + 1) / 2 + 1
            n_total = n_choose_k(n_col, k)
            n_limit = min(n_total, max_subsets)
            if (n_total > max_subsets) &
                print "('C(',i0,',',i0,') = ',i0,' > max_subsets = ',i0, &
                    &'; analyzing first ',i0,' subsets')", &
                    n_col, k, n_total, max_subsets, n_limit

            combo = [(j, j=1,k)]
            do n_sub = 1, n_limit
                print "(/,'=== subset (',i0,' assets): ',a)", k, &
                    trim(join(df_ret%columns(combo), " "))
                if (do_cov) then
                    print "(a)", "--- covariance changepoints ---"
                    call solve_changepoints(max_m, &
                        multivar_cost_matrix(R(:, combo), min_seg_len=min_seg_len), &
                        dp_table, parent)
                    call print_model_selection(dp_table, parent, best_bic_cp=best_bic_sub, &
                        params_per_seg=pps, print_each=print_all_models)
                    call corrmat_model_range(print_segs, best_bic_sub, max_m, m_lo, m_hi)
                    do ms = m_lo, m_hi
                        if (dp_table(n, ms) >= 1.0e19_dp) cycle
                        seg_ends(1:ms) = segment_ends(parent, ms)
                        call print_covmat_model(best_bic_sub, seg_ends(1:ms), &
                            R(:,combo), df_ret%columns(combo), ret_dates)
                    end do
                end if
                if (do_corr) then
                    print "(a)", "--- correlation changepoints ---"
                    call solve_changepoints(max_m, &
                        multivar_cost_matrix(R_std(:, combo), min_seg_len=min_seg_len), &
                        dp_table, parent)
                    call print_model_selection(dp_table, parent, best_bic_cp=best_bic_sub, &
                        params_per_seg=pps, print_each=print_all_models)
                    call corrmat_model_range(print_segs, best_bic_sub, max_m, m_lo, m_hi)
                    do ms = m_lo, m_hi
                        if (dp_table(n, ms) >= 1.0e19_dp) cycle
                        seg_ends(1:ms) = segment_ends(parent, ms)
                        call print_corrmat_model(best_bic_sub, seg_ends(1:ms), &
                            R(:,combo), df_ret%columns(combo), ret_dates, &
                            scale_ret, print_diffs, alpha_diff)
                    end do
                end if
                if (n_sub < n_limit) call next_combination(combo, n_col)
            end do
        end do
    end if

    deallocate(dp_table, parent, R, seg_ends)
    if (allocated(cp_aic)) deallocate(cp_aic, cp_bic)
    if (allocated(R_std))  deallocate(R_std)
    call print_wall_time(t_start)

end program xreturns_corrsub_cp
