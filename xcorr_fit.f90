!> Program: xcorr_fit
!!
!! Reads a bivariate time series from a file (sim_data.txt) and detects
!! changepoints in the correlation between the two series.
!!
!! The input file format should be:
!! # n
!! # x y rho
!! x1 y1 rho1
!! ...
!!
!! It uses dynamic programming to detect up to max_m changepoints,
!! outputting model selection statistics (AIC, BIC) and the
!! detected changepoint locations.

program xcorr_fit
    use kind_mod, only: dp, long_int
    use changepoint_mod, only: cost_matrix, solve_changepoints
    use util_mod, only: print_wall_time
    use io_utils_mod, only: print_model_selection, print_estimated_parameters
    implicit none

    integer, parameter :: max_m = 11, min_seg_len = 50
    integer :: n, iu, i
    real(kind=dp), allocatable :: x(:), y(:), cost(:,:), dp_table(:,:)
    integer, allocatable :: parent(:,:)
    character(len=10), allocatable :: date_labels(:)
    character(len=100) :: line
    integer(kind=long_int) :: t_start
    character (len=*), parameter :: data_file="sim_data.txt"
    call system_clock(t_start)
    
    ! Open file and read n
    open(newunit=iu, file=data_file, status="old")
    read(iu, "(a)") line ! Read "# n" line
    read(line(2:), *) n   ! Parse n
    read(iu, "(a)") line ! Skip header "# x y rho"

    allocate(x(n), y(n), cost(n, n), dp_table(n, max_m), parent(n, max_m), date_labels(n))

    do i = 1, n
        read(iu, *) x(i), y(i) ! rho is in the file but not needed for fitting
        write(date_labels(i), "(i0, '-01-01')") 2000 + (i-1)/365
    end do
    close(iu)
    print "('read ',i0,' observations from ',a)", n, trim(data_file)

    cost = cost_matrix(x, y, min_seg_len=min_seg_len)
    call solve_changepoints(max_m, cost, dp_table, parent)
    call print_model_selection(dp_table, parent)
    call print_estimated_parameters(max_m, parent, x, y)
    deallocate(x, y, cost, dp_table, parent, date_labels)
    call print_wall_time(t_start)
end program xcorr_fit
