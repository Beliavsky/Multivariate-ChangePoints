!> Program: xcorr_sim
!!
!! Simulates a bivariate normal time series with segments where 
!! the correlation parameter rho shifts.
!!
!! It then uses dynamic programming to detect changepoints in 
!! correlation by maximizing the log-likelihood of the segments.
!! Finally, it outputs model selection statistics (AIC, BIC) and 
!! the detected changepoint locations.

program xcorr_sim
    use kind_mod, only: dp, long_int
    use random_mod, only: rnorm
    use basic_stats_mod, only: mean
    use changepoint_mod, only: log_likelihood_corr, cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time, sort_int, set_segment_values
    use sim_changepoint_mod, only: generate_series
    use io_utils_mod, only: print_true_params, print_model_selection, print_estimated_parameters
    implicit none

    integer, parameter :: n = 3000, max_cp = 20, max_m = max_cp+1, nseg_true = 3, true_cps(nseg_true-1) = [1000, 2000]
    integer, parameter :: min_seg_len = 50
    real(kind=dp) :: x(n), y(n), corr_true(nseg_true) = [0.2_dp, 0.7_dp, -0.4_dp]
    real(kind=dp) :: rho_ts(n)
    real(kind=dp), allocatable :: dp_table(:,:), cost(:,:)
    integer, allocatable :: parent(:,:)
    integer(kind=long_int) :: t_start
    integer :: i, iu
    character (len=*), parameter :: sim_data_file = "sim_data.txt"
    logical, parameter :: write_sim = .true.
    character(len=10), allocatable :: date_labels(:)
    
    allocate(date_labels(n))
    do i = 1, n
        write(date_labels(i), "(i0, '-01-01')") 2000 + (i-1)/365
    end do
    
    call system_clock(t_start)
    allocate(dp_table(n, max_m), parent(n, max_m), cost(n, n))
    call generate_series(true_cps, corr_true, x, y)
    
    ! Determine true rho time series
    rho_ts = set_segment_values(n, true_cps, corr_true)

    if (write_sim) then
        ! Write simulated data
        open(newunit=iu, file=sim_data_file, status='replace')
        write(iu, "(a, i0)") "# ", n
        write(iu, "(a)") "# x y rho"
        do i = 1, n
            write(iu, "(3f18.12)") x(i), y(i), rho_ts(i)
        end do
        close(iu)
    end if

    call print_true_params(true_cps, corr_true, x, y)
    
    cost = cost_matrix(x, y, min_seg_len=min_seg_len)
    call solve_changepoints(max_m, cost, dp_table, parent)
    call print_model_selection(dp_table, parent)
    call print_estimated_parameters(max_m, parent, x, y)
    
    if (write_sim) print "(a)", "wrote simulated data to " // trim(sim_data_file)
    call print_wall_time(t_start)
    deallocate(cost, dp_table, parent, date_labels)
end program xcorr_sim
