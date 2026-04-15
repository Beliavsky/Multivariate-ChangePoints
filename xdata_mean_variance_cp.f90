!> Program: xdata_mean_variance_cp
!!
!! Reads a matrix data file and for each column finds changepoints
!! in the mean and/or variance using a profile normal mean-shift model.
!!
!! Input file format:
!!   First line:  # nrow ncol [any other text]
!!   Further lines beginning with '#' are treated as comments and skipped.
!!   Remaining lines: nrow rows of ncol whitespace-separated real values.
!!
!!   do_mean     = .true.  applies the cost to z_t = x_t        (mean shifts)
!!   do_variance = .true.  applies the cost to z_t = x_t^2      (if .not. use_log_var)
!!                         or z_t = log(var_offset + x_t^2)     (if use_log_var)
!!
!! print_segs: 0 = AIC/BIC table only
!!             1 = segment details for BIC-chosen model only
!!             2 = segment details for 0 changepoints through BIC-chosen
!!             3 = segment details for all models studied

program xdata_mean_variance_cp
    use kind_mod, only: dp, long_int
    use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints, segment_ends
    use util_mod, only: print_wall_time, read_matrix
    use io_utils_mod, only: print_model_selection, print_univar_segments, &
        corrmat_model_range
    implicit none

    integer, parameter :: max_cp = 100, max_m = max_cp + 1, min_seg_len = 20, max_cols = 10000
    character(len=*), parameter :: data_file = "sim_matrix_300_obs.txt"
    logical, parameter :: do_mean     = .false.
    logical, parameter :: do_variance = .true.
    integer, parameter :: print_segs  = 1       ! 0=none, 1=BIC, 2=0..BIC, 3=all
    logical, parameter :: print_table = .false.  ! .false. suppresses M/LL/AIC/BIC table
    integer, parameter :: max_col = 3
    logical, parameter :: use_log_var = .true.
    real(dp), parameter :: var_offset  = 0.01_dp  ! floor in log(var_offset + x^2)

    integer :: n, n_col, j, ba, bb, ms, m_lo, m_hi
    real(dp), allocatable :: x(:,:), dp_table(:,:), z(:)
    integer, allocatable :: parent(:,:)
    integer, allocatable :: mean_aic(:), mean_bic(:), var_aic(:), var_bic(:)
    character(len=20), allocatable :: col_labels(:)
    character(len=10), allocatable :: obs_labels(:)
    integer(kind=long_int) :: t_start

    call system_clock(t_start)

    call read_matrix(data_file, x, ncol_max = max_col)
    n = size(x,1)
    n_col = size(x,2)
    print "(a,i0,a,i0,a)", "read ", n, " rows x ", n_col, " cols from " // trim(data_file)

    allocate(dp_table(n, max_m), parent(n, max_m), z(n))
    allocate(mean_aic(n_col), mean_bic(n_col), var_aic(n_col), var_bic(n_col))
    mean_aic = 0;  mean_bic = 0
    var_aic  = 0;  var_bic  = 0

    ! Build column labels: "x1", "x2", ...
    allocate(col_labels(n_col))
    do j = 1, n_col
        write(col_labels(j), "(a,i0)") "x", j
    end do

    ! Build observation labels right-justified in 10-char fields
    allocate(obs_labels(n))
    do j = 1, n
        write(obs_labels(j), "(i10)") j
    end do

    ! ── mean changepoints ────────────────────────────────────────────────────────
    if (do_mean) then
        do j = 1, n_col
            print "(/,a)", "changes in mean of " // trim(col_labels(j))
            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(x(:,j), min_seg_len=min_seg_len), dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, &
                params_per_seg=3, print_each=print_table)
            mean_aic(j) = ba
            mean_bic(j) = bb
            if (print_segs > 0) then
                call corrmat_model_range(print_segs - 1, bb, max_m, m_lo, m_hi)
                do ms = m_lo, m_hi
                    if (dp_table(n, ms) >= 1.0e19_dp) cycle
                    call print_univar_segments(bb, segment_ends(parent, ms), x(:,j), &
                        trim(col_labels(j)), obs_labels)
                end do
            end if
        end do
    end if

    ! ── variance changepoints ────────────────────────────────────────────────────
    if (do_variance) then
        do j = 1, n_col
            if (use_log_var) then
                z = log(var_offset + x(:,j)**2)
            else
                z = x(:,j)**2
            end if
            print "(/,a)", "changes in " // var_label() // " of " // trim(col_labels(j))
            call solve_changepoints(max_m, &
                mean_shift_cost_matrix(z, min_seg_len=min_seg_len), dp_table, parent)
            call print_model_selection(dp_table, parent, best_aic_cp=ba, best_bic_cp=bb, &
                params_per_seg=3, print_each=print_table)
            var_aic(j) = ba
            var_bic(j) = bb
            if (print_segs > 0) then
                call corrmat_model_range(print_segs - 1, bb, max_m, m_lo, m_hi)
                do ms = m_lo, m_hi
                    if (dp_table(n, ms) >= 1.0e19_dp) cycle
                    call print_univar_segments(bb, segment_ends(parent, ms), x(:,j), &
                        trim(col_labels(j)), obs_labels)
                end do
            end if
        end do
    end if

    ! ── summary ──────────────────────────────────────────────────────────────────
    if (do_mean) then
        print "(/,a)", "Summary: # mean changepoints chosen by AIC and BIC"
        print "(a10,2a8)", "Series", "AIC", "BIC"
        do j = 1, n_col
            print "(a10,2i8)", trim(col_labels(j)), mean_aic(j), mean_bic(j)
        end do
    end if

    if (do_variance) then
        print "(/,'Summary: # ',a,' changepoints chosen by AIC and BIC')", var_label()
        print "(a10,2a8)", "Series", "AIC", "BIC"
        do j = 1, n_col
            print "(a10,2i8)", trim(col_labels(j)), var_aic(j), var_bic(j)
        end do
    end if

    deallocate(x, dp_table, parent, z)
    deallocate(mean_aic, mean_bic, var_aic, var_bic)
    deallocate(col_labels, obs_labels)
    call print_wall_time(t_start)

contains

    pure function var_label() result(label)
        character(len=20) :: label
        if (use_log_var) then
            label = "log-variance"
        else
            label = "variance"
        end if
    end function var_label

end program xdata_mean_variance_cp
