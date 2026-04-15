!> Program: xresample_cp
!!
!! Reads asset prices, computes returns, and for each asset finds changepoints
!! in the mean and/or variance using a profile normal mean-shift model applied
!! to the actual return series and to n_resample bootstrap-resampled copies.
!!
!! Bootstrap resampling (sampling rows with replacement) destroys any serial
!! structure, so it provides a null distribution for the number of changepoints
!! that would be found by chance.  Comparing the actual BIC changepoint count
!! to this distribution gives an empirical p-value.
!!
!!   do_mean     = .true.  detects mean shifts     (z_t = r_t)
!!   do_variance = .true.  detects variance shifts
!!                         use_log_var=.false.: z_t = r_t^2             (good locations, inverted null test)
!!                         use_log_var=.true.:  z_t = log(var_offset+r_t^2)  (approximately Gaussian, correct null test)
!!
!! Output:
!!   For the actual returns: the usual M / LL / AIC / BIC table and, if
!!   print_segs > 0, segment statistics for the BIC-chosen model.
!!
!!   For resampled runs: no per-run output; only running statistics.
!!
!!   Summary table (per series, per analysis type):
!!     Actual   -- BIC changepoint count for actual returns
!!     RMean    -- mean BIC changepoint count over resampled runs
!!     RMax     -- maximum BIC changepoint count over resampled runs
!!     p-value  -- fraction of resampled runs with BIC count >= Actual
!!
!! print_segs (actual data only): 0=none, 1=BIC only, 2=0..BIC, 3=all

program xresample_cp
    use kind_mod, only: dp, long_int
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow, ncol, operator(*)
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints, segment_ends
    use util_mod, only: print_wall_time
    use io_utils_mod, only: print_model_selection, print_return_segments, &
        corrmat_model_range, keep_obs
    implicit none

    type(DataFrame_index_date) :: df_px, df_ret, df_shuf
    integer, parameter :: max_cp = 100, max_m = max_cp + 1, min_seg_len = 50, max_assets = 1000
    character(len=*), parameter :: prices_file = "asset_class_etf_prices.csv", &
        sym_allowed(*) = [character(len=5) :: "SPY"]
    real(kind=dp), parameter :: scale_ret = 100.0_dp
    integer, parameter :: max_days   = 0       ! 0 = use all returns
    logical, parameter :: latest     = .true.  ! if max_days > 0: .true. = latest
    logical, parameter :: do_mean     = .false.
    logical, parameter :: do_variance = .true.
    integer, parameter :: print_segs  = 1      ! for actual returns only: 0=none, 1=BIC, 2=0..BIC, 3=all
    integer, parameter :: n_resample  = 100    ! number of bootstrap resample runs
    ! Variance series: use_log_var=.true.  => z = log(var_offset + r^2), approximately Gaussian,
    !                                         corrects the inverted null-hypothesis test
    !                  use_log_var=.false. => z = r^2, finds good changepoint locations
    !                                         but resampled data finds more CPs than actual (test inverted)
    logical,       parameter :: use_log_var = .true.
    real(kind=dp), parameter :: var_offset  = 0.01_dp  ! floor in log(var_offset + r^2); ignored if .not. use_log_var

    integer :: n, n_col, j, ba, bb, ms, m_lo, m_hi, ir
    real(kind=dp), allocatable :: dp_table(:,:), R(:,:), z(:)
    integer, allocatable :: parent(:,:)
    integer, allocatable :: mean_bic_act(:), var_bic_act(:)
    integer, allocatable :: mean_bic_rsum(:), mean_bic_rmax(:), mean_bic_rcnt(:)
    integer, allocatable :: var_bic_rsum(:),  var_bic_rmax(:),  var_bic_rcnt(:)
    integer(kind=long_int) :: t_start
    character(len=10), allocatable :: ret_dates(:)

    call system_clock(t_start)

    call df_px%read_csv(prices_file, max_col=max_assets)
    if (size(sym_allowed) > 0) df_px = df_px%select(columns=sym_allowed)
    print "(*(1x,a,1x,i0))", "read", nrow(df_px), "days and", ncol(df_px), &
        "columns from " // trim(prices_file)
    if (max_days > 0) call keep_obs(df_px, max_days, latest, verbose=.true.)
    df_ret = scale_ret * df_px%pct_change(dropna=.true.)
    if (scale_ret /= 1.0_dp) print "('return scaling: ', f0.4)", scale_ret

    n     = size(df_ret%values, 1)
    n_col = ncol(df_ret)
    ret_dates = df_ret%index%to_str()

    allocate(R(n, n_col), dp_table(n, max_m), parent(n, max_m), z(n))
    R = df_ret%values

    allocate(mean_bic_act(n_col), var_bic_act(n_col))
    allocate(mean_bic_rsum(n_col), mean_bic_rmax(n_col), mean_bic_rcnt(n_col))
    allocate(var_bic_rsum(n_col),  var_bic_rmax(n_col),  var_bic_rcnt(n_col))
    mean_bic_act  = 0;  var_bic_act  = 0
    mean_bic_rsum = 0;  mean_bic_rmax = 0;  mean_bic_rcnt = 0
    var_bic_rsum  = 0;  var_bic_rmax  = 0;  var_bic_rcnt  = 0

    ! ── actual returns ──────────────────────────────────────────────────────────────
    print "(/,a)", "=== Actual returns ==="

    if (do_mean) then
        do j = 1, n_col
            print "(/,a)", "changes in mean of " // trim(df_ret%columns(j))
            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(R(:,j), min_seg_len=min_seg_len), dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, &
                params_per_seg=3)
            mean_bic_act(j) = bb
            if (print_segs > 0) then
                call corrmat_model_range(print_segs - 1, bb, max_m, m_lo, m_hi)
                do ms = m_lo, m_hi
                    if (dp_table(n, ms) >= 1.0e19_dp) cycle
                    call print_return_segments(bb, segment_ends(parent, ms), R(:,j), &
                        trim(df_ret%columns(j)), ret_dates, scale_ret)
                end do
            end if
        end do
    end if

    if (do_variance) then
        do j = 1, n_col
            if (use_log_var) then
                z = log(var_offset + R(:,j)**2)
            else
                z = R(:,j)**2
            end if
            print "(/,a)", "changes in " // var_label() // " of " // trim(df_ret%columns(j))
            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(z, min_seg_len=min_seg_len), dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, &
                params_per_seg=3)
            var_bic_act(j) = bb
            if (print_segs > 0) then
                call corrmat_model_range(print_segs - 1, bb, max_m, m_lo, m_hi)
                do ms = m_lo, m_hi
                    if (dp_table(n, ms) >= 1.0e19_dp) cycle
                    call print_return_segments(bb, segment_ends(parent, ms), R(:,j), &
                        trim(df_ret%columns(j)), ret_dates, scale_ret)
                end do
            end if
        end do
    end if

    ! ── resampled returns ───────────────────────────────────────────────────────────
    print "(/,'=== Resampling: ',i0,' runs ===')", n_resample

    do ir = 1, n_resample
        write(*, "(i0,' ')", advance='no') ir
        df_shuf = df_ret%resample()

        if (do_mean) then
            do j = 1, n_col
                call solve_changepoints(max_m, &
                    mean_shift_cost_matrix(df_shuf%values(:,j), min_seg_len=min_seg_len), &
                    dp_table, parent)
                bb = bic_best(dp_table, n, pps=3)
                mean_bic_rsum(j) = mean_bic_rsum(j) + bb
                if (bb > mean_bic_rmax(j)) mean_bic_rmax(j) = bb
                if (bb >= mean_bic_act(j)) mean_bic_rcnt(j) = mean_bic_rcnt(j) + 1
            end do
        end if

        if (do_variance) then
            do j = 1, n_col
                if (use_log_var) then
                    z = log(var_offset + df_shuf%values(:,j)**2)
                else
                    z = df_shuf%values(:,j)**2
                end if
                call solve_changepoints(max_m, &
                    mean_shift_cost_matrix(z, min_seg_len=min_seg_len), dp_table, parent)
                bb = bic_best(dp_table, n, pps=3)
                var_bic_rsum(j) = var_bic_rsum(j) + bb
                if (bb > var_bic_rmax(j)) var_bic_rmax(j) = bb
                if (bb >= var_bic_act(j)) var_bic_rcnt(j) = var_bic_rcnt(j) + 1
            end do
        end if
    end do

    ! ── summary ────────────────────────────────────────────────────────────────────
    if (do_mean) then
        print "(/,'Summary: mean changepoints (BIC), actual vs.',i0,' resampled runs')", n_resample
        print "(a10, a7, a9, a7, a9)", "Series", "Actual", "RMean", "RMax", "p-value"
        do j = 1, n_col
            print "(a10, i7, f9.2, i7, f9.3)", &
                trim(df_ret%columns(j)), mean_bic_act(j), &
                real(mean_bic_rsum(j), dp) / n_resample, &
                mean_bic_rmax(j), &
                real(mean_bic_rcnt(j), dp) / n_resample
        end do
    end if

    if (do_variance) then
        print "(/,'Summary: ',a,' changepoints (BIC), actual vs.',i0,' resampled runs')", &
            var_label(), n_resample
        print "(a10, a7, a9, a7, a9)", "Series", "Actual", "RMean", "RMax", "p-value"
        do j = 1, n_col
            print "(a10, i7, f9.2, i7, f9.3)", &
                trim(df_ret%columns(j)), var_bic_act(j), &
                real(var_bic_rsum(j), dp) / n_resample, &
                var_bic_rmax(j), &
                real(var_bic_rcnt(j), dp) / n_resample
        end do
    end if

    deallocate(R, dp_table, parent, z)
    deallocate(mean_bic_act, var_bic_act)
    deallocate(mean_bic_rsum, mean_bic_rmax, mean_bic_rcnt)
    deallocate(var_bic_rsum,  var_bic_rmax,  var_bic_rcnt)
    call print_wall_time(t_start)

contains

    pure function var_label() result(label)
        !> Human-readable name for the variance series currently in use.
        character(len=20) :: label
        if (use_log_var) then
            label = "log-variance"
        else
            label = "variance"
        end if
    end function var_label

    pure function bic_best(dp_tab, nt, pps) result(best_cp)
        !> BIC-optimal number of changepoints, computed directly from dp_table.
        real(kind=dp), intent(in) :: dp_tab(:,:)
        integer,       intent(in) :: nt, pps
        integer :: best_cp, m
        real(kind=dp) :: bic, min_bic
        min_bic = huge(1.0_dp)
        best_cp = 0
        do m = 1, size(dp_tab, 2)
            if (dp_tab(nt, m) >= 1.0e19_dp) cycle
            bic = 2.0_dp * dp_tab(nt, m) + real(pps * m - 1, dp) * log(real(nt, dp))
            if (bic < min_bic) then
                min_bic = bic
                best_cp = m - 1
            end if
        end do
    end function bic_best

end program xresample_cp
