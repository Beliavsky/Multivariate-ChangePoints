module sim_changepoint_mod
    use kind_mod, only: dp
    use random_mod, only: rnorm
    implicit none
    public :: generate_series, generate_multivar_series

contains

    subroutine generate_series(true_cps, corr_true, x, y)
        !> Simulates a bivariate normal time series with shifts in correlation.
        integer, intent(in) :: true_cps(:)        ! Changepoint locations
        real(kind=dp), intent(in) :: corr_true(:) ! Correlation coefficients per segment
        real(kind=dp), intent(out) :: x(:), y(:)  ! Output series (n)
        integer :: i, k, n, nseg_true

        if (size(x) /= size(y)) &
            error stop "generate_series: size(x) /= size(y)"
        if (size(corr_true) /= size(true_cps) + 1) &
            error stop "generate_series: size(corr_true) /= size(true_cps) + 1"
        
        n = size(x)
        nseg_true = size(corr_true)
        x = rnorm(n, 0.0_dp, 1.0_dp)
        y = rnorm(n, 0.0_dp, 1.0_dp)
        do i = 1, n
            k = nseg_true
            do k = 1, nseg_true - 1
                if (i <= true_cps(k)) then
                    exit
                end if
            end do
            y(i) = corr_true(k)*x(i) + sqrt(max(0.0_dp, 1.0_dp-corr_true(k)**2))*y(i)
        end do
    end subroutine generate_series

    !> Cholesky factorisation: A = L * L^T (lower triangle only).
    !! ok = .false. if A is not positive definite.
    subroutine cholesky_lower(A, L, ok)
        real(kind=dp), intent(in)  :: A(:,:)
        real(kind=dp), intent(out) :: L(:,:)
        logical,       intent(out) :: ok
        integer :: p, i, j
        real(kind=dp) :: s
        p = size(A, 1)
        L = 0.0_dp
        ok = .true.
        do j = 1, p
            s = A(j,j) - sum(L(j, 1:j-1)**2)
            if (s <= 0.0_dp) then
                ok = .false.
                return
            end if
            L(j,j) = sqrt(s)
            do i = j+1, p
                L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
            end do
        end do
    end subroutine cholesky_lower

    !> Simulate an n×p multivariate normal series with piecewise-constant covariance.
    !! cov_true(:,:,k) is the p×p covariance matrix for segment k.
    !! true_cps(k) is the last observation of segment k (1..nseg-1).
    subroutine generate_multivar_series(true_cps, cov_true, R)
        integer,       intent(in)  :: true_cps(:)    ! nseg-1 breakpoints
        real(kind=dp), intent(in)  :: cov_true(:,:,:)! p×p×nseg
        real(kind=dp), intent(out) :: R(:,:)         ! n×p
        integer :: n, p, nseg, i, k, seg, seg_prev
        real(kind=dp) :: z(size(R,2)), L(size(R,2), size(R,2))
        logical :: ok

        n    = size(R, 1)
        p    = size(R, 2)
        nseg = size(cov_true, 3)

        seg_prev = 0
        do i = 1, n
            seg = nseg
            do k = 1, nseg - 1
                if (i <= true_cps(k)) then
                    seg = k
                    exit
                end if
            end do
            if (seg /= seg_prev) then          ! recompute Cholesky at segment changes
                call cholesky_lower(cov_true(:,:,seg), L, ok)
                if (.not. ok) error stop "generate_multivar_series: covariance not positive definite"
                seg_prev = seg
            end if
            z = rnorm(p, 0.0_dp, 1.0_dp)
            R(i,:) = matmul(L, z)
        end do
    end subroutine generate_multivar_series

end module sim_changepoint_mod
