module changepoint_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: log_likelihood_corr, cost_matrix, cost_matrix_slow, mean_shift_cost_matrix, &
              multivar_cost_matrix, solve_changepoints
    public :: solve_continuous_pwl_sse, fit_continuous_pwl_given_cps

    type :: quad_state
        integer :: nquad = 0
        real(kind=dp), allocatable :: a(:), b(:), c(:)
        integer, allocatable :: parent_t(:), parent_q(:)
    end type quad_state

contains

    !> Computes the negative log-likelihood of bivariate normal data segment.
    pure function log_likelihood_corr(x, y, min_len) result(ll)
        real(kind=dp), intent(in) :: x(:), y(:)
        integer, intent(in), optional :: min_len
        real(kind=dp) :: ll, r, ma, mb, va, vb
        integer :: n, min_len_

        min_len_ = 50
        if (present(min_len)) min_len_ = min_len

        n = size(x)
        if (n < min_len_) then
            ll = -1.0e20_dp
            return
        end if
        ma = sum(x)/n
        mb = sum(y)/n
        va = sum((x-ma)**2) / n
        vb = sum((y-mb)**2) / n
        r = sum((x-ma)*(y-mb)) / (sqrt(max(1e-20_dp, va * vb)) * n)

        ll = -0.5_dp * n * log(max(1.0e-10_dp, 1.0_dp - min(r**2, 1.0_dp - 1.0e-10_dp)))
    end function log_likelihood_corr

    !> Constructs a cost matrix for all possible segments using negative log-likelihood (O(N^2)).
    function cost_matrix(x, y, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: x(:), y(:)
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(x), size(x))
        real(kind=dp) :: sx(size(x)), sy(size(x)), sxx(size(x)), syy(size(x)), sxy(size(x))
        integer :: n, i, j, min_len

        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        n = size(x)

        sx(1) = x(1); sy(1) = y(1); sxx(1) = x(1)**2; syy(1) = y(1)**2; sxy(1) = x(1)*y(1)
        do i = 2, n
            sx(i) = sx(i-1) + x(i)
            sy(i) = sy(i-1) + y(i)
            sxx(i) = sxx(i-1) + x(i)**2
            syy(i) = syy(i-1) + y(i)**2
            sxy(i) = sxy(i-1) + x(i)*y(i)
        end do

        do i = 1, n
            do j = i, n
                if (j - i + 1 < min_len) then
                    cost(i, j) = 1.0e20_dp
                else
                    call fast_ll(i, j, sx, sy, sxx, syy, sxy, cost(i, j))
                end if
            end do
        end do
        do i = 2, n
            do j = 1, i - 1
                cost(i, j) = 1.0e20_dp
            end do
        end do
    end function cost_matrix

    subroutine fast_ll(i, j, sx, sy, sxx, syy, sxy, cost_val)
        integer, intent(in) :: i, j
        real(kind=dp), intent(in) :: sx(:), sy(:), sxx(:), syy(:), sxy(:)
        real(kind=dp), intent(out) :: cost_val
        real(kind=dp) :: n, r, ma, mb, va, vb, sx_s, sy_s, sxx_s, syy_s, sxy_s

        n = real(j - i + 1, dp)
        if (i == 1) then
            sx_s = sx(j); sy_s = sy(j); sxx_s = sxx(j); syy_s = syy(j); sxy_s = sxy(j)
        else
            sx_s = sx(j) - sx(i-1)
            sy_s = sy(j) - sy(i-1)
            sxx_s = sxx(j) - sxx(i-1)
            syy_s = syy(j) - syy(i-1)
            sxy_s = sxy(j) - sxy(i-1)
        end if
        ma = sx_s / n
        mb = sy_s / n
        va = (sxx_s - n * ma**2) / n
        vb = (syy_s - n * mb**2) / n
        r = (sxy_s - n * ma * mb) / (sqrt(max(1e-20_dp, va * vb)) * n)

        cost_val = 0.5_dp * n * log(max(1.0e-10_dp, 1.0_dp - min(r**2, 1.0_dp - 1.0e-10_dp)))
    end subroutine fast_ll

    !> Cost matrix for a 1-D series z using the profile normal mean-shift log-likelihood.
    !! cost(i,j) = (j-i+1)/2 * log(max(sample_variance(z(i:j)), eps))
    !! Use z = x*y for covariance changepoints; z = x*x for variance changepoints.
    function mean_shift_cost_matrix(z, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(z), size(z))
        real(kind=dp) :: sz(size(z)), szz(size(z))
        real(kind=dp) :: m, sz_s, szz_s, s2
        integer :: n, i, j, min_len

        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        n = size(z)

        sz(1) = z(1);  szz(1) = z(1)**2
        do i = 2, n
            sz(i)  = sz(i-1)  + z(i)
            szz(i) = szz(i-1) + z(i)**2
        end do

        do i = 1, n
            do j = i, n
                if (j - i + 1 < min_len) then
                    cost(i, j) = 1.0e20_dp
                else
                    m     = real(j - i + 1, dp)
                    sz_s  = sz(j)  - merge(sz(i-1),  0.0_dp, i > 1)
                    szz_s = szz(j) - merge(szz(i-1), 0.0_dp, i > 1)
                    s2    = szz_s / m - (sz_s / m)**2
                    cost(i, j) = 0.5_dp * m * log(max(s2, 1.0e-20_dp))
                end if
            end do
            do j = 1, i - 1
                cost(i, j) = 1.0e20_dp
            end do
        end do
    end function mean_shift_cost_matrix

    !> Log-determinant of a symmetric positive-definite matrix via Cholesky.
    !! Returns -1e30 if the matrix is not positive definite.
    pure function log_det_chol(A) result(ld)
        real(kind=dp), intent(in) :: A(:,:)
        real(kind=dp) :: ld
        integer :: p, i, j
        real(kind=dp) :: s, L(size(A,1), size(A,1))
        p = size(A, 1)
        L = 0.0_dp
        do j = 1, p
            s = A(j,j) - sum(L(j, 1:j-1)**2)
            if (s <= 0.0_dp) then
                ld = -1.0e30_dp
                return
            end if
            L(j,j) = sqrt(s)
            do i = j+1, p
                L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
            end do
        end do
        ld = 2.0_dp * sum([(log(L(i,i)), i=1,p)])
    end function log_det_chol

    !> Joint covariance-matrix changepoint cost matrix.
    !! cost(i,j) = m/2 * log|Sigma_hat(i:j)|  where Sigma_hat is the biased sample
    !! covariance matrix of the p-column return matrix R over rows i..j.
    !! Uses O(n*p^2) prefix sums; each cell costs O(p^2) to assemble + O(p^3) Cholesky.
    function multivar_cost_matrix(R, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: R(:,:)      ! n × p
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(R,1), size(R,1))
        real(kind=dp), allocatable :: SR(:,:), SP(:,:,:)
        real(kind=dp) :: S(size(R,2), size(R,2)), sr_s(size(R,2)), m_r, ld
        integer :: n, p, i, j, a, b, min_len

        n = size(R, 1);  p = size(R, 2)
        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        allocate(SR(0:n, p), SP(0:n, p, p))
        cost = 1.0e20_dp

        SR(0,:) = 0.0_dp;  SP(0,:,:) = 0.0_dp
        do i = 1, n
            SR(i,:) = SR(i-1,:) + R(i,:)
            do a = 1, p
                SP(i,a,:) = SP(i-1,a,:) + R(i,a) * R(i,:)
            end do
        end do

        do j = min_len, n
            do i = 1, j - min_len + 1
                m_r  = real(j - i + 1, dp)
                sr_s = SR(j,:) - SR(i-1,:)
                do a = 1, p
                    do b = a, p
                        S(a,b) = (SP(j,a,b) - SP(i-1,a,b)) / m_r &
                               - (sr_s(a)/m_r) * (sr_s(b)/m_r)
                        S(b,a) = S(a,b)
                    end do
                end do
                ld = log_det_chol(S)
                if (ld > -1.0e19_dp) cost(i,j) = 0.5_dp * m_r * ld
            end do
        end do
    end function multivar_cost_matrix

    !> Constructs a cost matrix for all possible segments using negative log-likelihood (slow O(N^3) version).
    pure function cost_matrix_slow(x, y) result(cost)
        real(kind=dp), intent(in) :: x(:), y(:)
        real(kind=dp) :: cost(size(x), size(x))
        integer :: n, i, j
        n = size(x)
        do i = 1, n
            do j = i, n
                cost(i, j) = -log_likelihood_corr(x(i:j), y(i:j))
            end do
        end do
    end function cost_matrix_slow

    subroutine solve_changepoints(max_m, cost, dp_table, parent)
        !> Solves the changepoint problem using dynamic programming to minimize total cost.
        integer, intent(in) :: max_m
        real(kind=dp), intent(in) :: cost(:, :)
        real(kind=dp), intent(out) :: dp_table(size(cost, 1), max_m)
        integer, intent(out) :: parent(size(cost, 1), max_m)
        integer :: i, m, k, n

        n = size(cost, 1)

        dp_table = 1.0e20_dp
        parent = 0
        do i = 1, n
            dp_table(i, 1) = cost(1, i)
        end do

        do m = 2, max_m
            do i = m, n
                do k = m-1, i-1
                    if (dp_table(k, m-1) + cost(k+1, i) < dp_table(i, m)) then
                        dp_table(i, m) = dp_table(k, m-1) + cost(k+1, i)
                        parent(i, m) = k
                    end if
                end do
            end do
        end do
    end subroutine solve_changepoints

    !> Solve the exact continuous piecewise linear least-squares problem for y.
    !!
    !! The fitted curve is continuous and piecewise linear, with observations 1:n.
    !! For m segments there are m-1 internal changepoints. Segment 1 fits y(1:t1).
    !! Later segments fit y(t_prev+1:t_curr), so each observation is used exactly once.
    !!
    !! The dynamic programming state for each (m,t) is represented as a lower envelope
    !! of quadratics in the endpoint value at time t. This gives an exact optimizer for
    !! the least-squares continuous PWL objective, not the greedy knot-insertion heuristic.
    subroutine solve_continuous_pwl_sse(y, max_m, min_seg_len, sse, cp_store, n_models_fitted)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: max_m
        integer, intent(in), optional :: min_seg_len
        real(kind=dp), intent(out) :: sse(max_m)
        integer, intent(out) :: cp_store(max(1, max_m - 1), max_m)
        integer, intent(out) :: n_models_fitted

        type(quad_state), allocatable :: states(:, :)
        real(kind=dp), allocatable :: sy(:), siy(:), syy(:)
        real(kind=dp), allocatable :: raw_a(:), raw_b(:), raw_c(:)
        integer, allocatable :: raw_pt(:), raw_pq(:), keep(:)
        integer :: n, min_len, max_m_eff, m, t, s, j, nraw, nkeep, best_q
        real(kind=dp) :: a0, b0, c0, qa, qb, qc, qd, qe, qf
        real(kind=dp) :: best_cost, denom, lin

        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        n = size(y)
        if (max_m < 1) error stop 'solve_continuous_pwl_sse: max_m must be positive'
        if (n < min_len) error stop 'solve_continuous_pwl_sse: size(y) < min_seg_len'

        max_m_eff = min(max_m, n / min_len)
        sse = huge(1.0_dp)
        cp_store = 0
        n_models_fitted = 0

        call build_prefix_sums(y, sy, siy, syy)
        allocate(states(max_m_eff, n))

        do t = min_len, n
            call first_segment_quadratic(t, sy, siy, syy, qa, qb, qc)
            call set_state_single(states(1, t), qa, qb, qc)
        end do

        do m = 2, max_m_eff
            do t = m * min_len, n
                nraw = 0
                do s = (m - 1) * min_len, t - min_len
                    nraw = nraw + states(m - 1, s)%nquad
                end do
                if (nraw == 0) cycle

                allocate(raw_a(nraw), raw_b(nraw), raw_c(nraw), raw_pt(nraw), raw_pq(nraw))
                nraw = 0
                do s = (m - 1) * min_len, t - min_len
                    if (states(m - 1, s)%nquad == 0) cycle
                    call later_segment_coefficients(s, t, sy, siy, syy, qa, qb, qc, qd, qe, qf)
                    do j = 1, states(m - 1, s)%nquad
                        a0 = states(m - 1, s)%a(j)
                        b0 = states(m - 1, s)%b(j)
                        c0 = states(m - 1, s)%c(j)
                        denom = a0 + qa
                        lin = b0 + qd
                        nraw = nraw + 1
                        raw_a(nraw) = qc - (qb * qb) / denom
                        raw_b(nraw) = qe - qb * lin / denom
                        raw_c(nraw) = c0 + qf - 0.25_dp * lin * lin / denom
                        raw_pt(nraw) = s
                        raw_pq(nraw) = j
                    end do
                end do

                call keep_lower_envelope(raw_a, raw_b, raw_c, keep, nkeep)
                call set_state_subset(states(m, t), raw_a, raw_b, raw_c, raw_pt, raw_pq, keep, nkeep)
                deallocate(raw_a, raw_b, raw_c, raw_pt, raw_pq, keep)
            end do
        end do

        do m = 1, max_m_eff
            if (states(m, n)%nquad == 0) cycle
            call best_terminal_quadratic(states(m, n)%a, states(m, n)%b, states(m, n)%c, best_q, best_cost)
            sse(m) = best_cost
            call backtrack_pwl_cps(states, m, n, best_q, cp_store(1:m-1, m))
            n_models_fitted = m
        end do

        deallocate(states, sy, siy, syy)
    end subroutine solve_continuous_pwl_sse

    !> Fit a continuous piecewise linear regression to y for a fixed changepoint set.
    subroutine fit_continuous_pwl_given_cps(y, cps, fitted, rss)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(out) :: fitted(:)
        real(kind=dp), intent(out), optional :: rss

        integer :: n, p, i, k
        real(kind=dp), allocatable :: xtx(:, :), xty(:), beta(:)
        real(kind=dp) :: t

        n = size(y)
        if (size(fitted) /= n) error stop 'fit_continuous_pwl_given_cps: size(fitted) /= size(y)'

        p = size(cps) + 2
        allocate(xtx(p, p), xty(p), beta(p))
        xtx = 0.0_dp
        xty = 0.0_dp

        do i = 1, n
            t = real(i, dp)
            call accumulate_normal_equations(t, y(i), cps, xtx, xty)
        end do

        call solve_linear_system(xtx, xty, beta)

        do i = 1, n
            t = real(i, dp)
            fitted(i) = beta(1) + beta(2) * t
            do k = 1, size(cps)
                fitted(i) = fitted(i) + beta(k + 2) * max(0.0_dp, t - real(cps(k), dp))
            end do
        end do

        if (present(rss)) rss = sum((y - fitted)**2)

        deallocate(xtx, xty, beta)
    end subroutine fit_continuous_pwl_given_cps

    subroutine build_prefix_sums(y, sy, siy, syy)
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), allocatable, intent(out) :: sy(:), siy(:), syy(:)
        integer :: n, i

        n = size(y)
        allocate(sy(0:n), siy(0:n), syy(0:n))
        sy = 0.0_dp
        siy = 0.0_dp
        syy = 0.0_dp
        do i = 1, n
            sy(i) = sy(i-1) + y(i)
            siy(i) = siy(i-1) + real(i, dp) * y(i)
            syy(i) = syy(i-1) + y(i) * y(i)
        end do
    end subroutine build_prefix_sums

    subroutine first_segment_quadratic(t, sy, siy, syy, a, b, c)
        integer, intent(in) :: t
        real(kind=dp), intent(in) :: sy(0:), siy(0:), syy(0:)
        real(kind=dp), intent(out) :: a, b, c
        real(kind=dp) :: qa, qb, qc, qd, qe, qf, d, sumy, sumiy, sumjy, ya, yb, denom

        d = real(t - 1, dp)
        sumy = sy(t)
        sumiy = siy(t)
        sumjy = sumiy - sumy

        qa = (d + 1.0_dp) * (2.0_dp * d + 1.0_dp) / (6.0_dp * d)
        qb = (d * d - 1.0_dp) / (6.0_dp * d)
        qc = qa
        ya = sumy - sumjy / d
        yb = sumjy / d
        qd = -2.0_dp * ya
        qe = -2.0_dp * yb
        qf = syy(t)

        denom = qa
        a = qc - (qb * qb) / denom
        b = qe - qb * qd / denom
        c = qf - 0.25_dp * qd * qd / denom
    end subroutine first_segment_quadratic

    subroutine later_segment_coefficients(s, t, sy, siy, syy, a, b, c, dcoef, ecoef, fcoef)
        integer, intent(in) :: s, t
        real(kind=dp), intent(in) :: sy(0:), siy(0:), syy(0:)
        real(kind=dp), intent(out) :: a, b, c, dcoef, ecoef, fcoef
        real(kind=dp) :: len, sumy, sumiy, sumjy, ya, yb

        len = real(t - s, dp)
        sumy = sy(t) - sy(s)
        sumiy = siy(t) - siy(s)
        sumjy = sumiy - real(s, dp) * sumy

        a = (len - 1.0_dp) * (2.0_dp * len - 1.0_dp) / (6.0_dp * len)
        b = (len * len - 1.0_dp) / (6.0_dp * len)
        c = (len + 1.0_dp) * (2.0_dp * len + 1.0_dp) / (6.0_dp * len)
        ya = sumy - sumjy / len
        yb = sumjy / len
        dcoef = -2.0_dp * ya
        ecoef = -2.0_dp * yb
        fcoef = syy(t) - syy(s)
    end subroutine later_segment_coefficients

    subroutine set_state_single(state, a, b, c)
        type(quad_state), intent(inout) :: state
        real(kind=dp), intent(in) :: a, b, c

        state%nquad = 1
        allocate(state%a(1), state%b(1), state%c(1), state%parent_t(1), state%parent_q(1))
        state%a(1) = a
        state%b(1) = b
        state%c(1) = c
        state%parent_t(1) = 0
        state%parent_q(1) = 0
    end subroutine set_state_single

    subroutine set_state_subset(state, a, b, c, parent_t, parent_q, keep, nkeep)
        type(quad_state), intent(inout) :: state
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        integer, intent(in) :: parent_t(:), parent_q(:), keep(:), nkeep
        integer :: i, idx

        state%nquad = nkeep
        allocate(state%a(nkeep), state%b(nkeep), state%c(nkeep), state%parent_t(nkeep), state%parent_q(nkeep))
        do i = 1, nkeep
            idx = keep(i)
            state%a(i) = a(idx)
            state%b(i) = b(idx)
            state%c(i) = c(idx)
            state%parent_t(i) = parent_t(idx)
            state%parent_q(i) = parent_q(idx)
        end do
    end subroutine set_state_subset

    subroutine best_terminal_quadratic(a, b, c, best_q, best_cost)
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        integer, intent(out) :: best_q
        real(kind=dp), intent(out) :: best_cost
        integer :: j
        real(kind=dp) :: val

        best_q = 1
        best_cost = c(1) - 0.25_dp * b(1) * b(1) / a(1)
        do j = 2, size(a)
            val = c(j) - 0.25_dp * b(j) * b(j) / a(j)
            if (val < best_cost) then
                best_cost = val
                best_q = j
            end if
        end do
    end subroutine best_terminal_quadratic

    subroutine backtrack_pwl_cps(states, m, t, q, cps)
        type(quad_state), intent(in) :: states(:, :)
        integer, intent(in) :: m, t, q
        integer, intent(out) :: cps(:)
        integer :: seg, cur_t, cur_q

        if (size(cps) /= max(0, m - 1)) error stop 'backtrack_pwl_cps: invalid cps size'

        cur_t = t
        cur_q = q
        do seg = m, 2, -1
            cps(seg - 1) = states(seg, cur_t)%parent_t(cur_q)
            cur_q = states(seg, cur_t)%parent_q(cur_q)
            cur_t = cps(seg - 1)
        end do
    end subroutine backtrack_pwl_cps

    subroutine keep_lower_envelope(a, b, c, keep, nkeep)
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        integer, allocatable, intent(out) :: keep(:)
        integer, intent(out) :: nkeep

        integer :: k, i, j, nroot, nuniq, idx, nsurvive
        real(kind=dp), allocatable :: roots(:), roots_u(:), aa_s(:), bb_s(:), cc_s(:)
        logical, allocatable :: active(:), survive(:)
        real(kind=dp) :: aa, bb, cc, disc, r1, r2, x, delta

        k = size(a)
        if (k == 0) then
            allocate(keep(0))
            nkeep = 0
            return
        end if

        allocate(survive(k))
        survive = .true.
        do i = 1, k
            do j = 1, k
                if (i == j) cycle
                if (quadratic_dominated(a(i), b(i), c(i), a(j), b(j), c(j))) then
                    survive(i) = .false.
                    exit
                end if
            end do
        end do

        nsurvive = count(survive)
        allocate(aa_s(nsurvive), bb_s(nsurvive), cc_s(nsurvive), keep(nsurvive))
        idx = 0
        do i = 1, k
            if (survive(i)) then
                idx = idx + 1
                aa_s(idx) = a(i)
                bb_s(idx) = b(i)
                cc_s(idx) = c(i)
                keep(idx) = i
            end if
        end do
        deallocate(survive)

        allocate(active(nsurvive))
        active = .false.
        allocate(roots(max(1, nsurvive * (nsurvive - 1))))
        nroot = 0

        do i = 1, nsurvive - 1
            do j = i + 1, nsurvive
                aa = aa_s(i) - aa_s(j)
                bb = bb_s(i) - bb_s(j)
                cc = cc_s(i) - cc_s(j)
                if (abs(aa) <= 1.0e-12_dp) then
                    if (abs(bb) > 1.0e-12_dp) then
                        nroot = nroot + 1
                        roots(nroot) = -cc / bb
                    end if
                else
                    disc = bb * bb - 4.0_dp * aa * cc
                    if (disc >= -1.0e-10_dp) then
                        disc = max(disc, 0.0_dp)
                        r1 = (-bb - sqrt(disc)) / (2.0_dp * aa)
                        r2 = (-bb + sqrt(disc)) / (2.0_dp * aa)
                        nroot = nroot + 1
                        roots(nroot) = r1
                        if (abs(r2 - r1) > 1.0e-10_dp) then
                            nroot = nroot + 1
                            roots(nroot) = r2
                        end if
                    end if
                end if
            end do
        end do

        if (nroot > 0) then
            call sort_real(roots(1:nroot))
            allocate(roots_u(nroot))
            nuniq = 1
            roots_u(1) = roots(1)
            do i = 2, nroot
                if (abs(roots(i) - roots_u(nuniq)) > 1.0e-8_dp) then
                    nuniq = nuniq + 1
                    roots_u(nuniq) = roots(i)
                end if
            end do

            do i = 1, nuniq
                call mark_minimizers(roots_u(i), aa_s, bb_s, cc_s, active)
            end do

            delta = abs(roots_u(1)) + 1.0_dp
            call mark_minimizers(roots_u(1) - delta, aa_s, bb_s, cc_s, active)
            do i = 1, nuniq - 1
                x = 0.5_dp * (roots_u(i) + roots_u(i + 1))
                call mark_minimizers(x, aa_s, bb_s, cc_s, active)
            end do
            delta = abs(roots_u(nuniq)) + 1.0_dp
            call mark_minimizers(roots_u(nuniq) + delta, aa_s, bb_s, cc_s, active)
            deallocate(roots_u)
        else
            call mark_minimizers(0.0_dp, aa_s, bb_s, cc_s, active)
        end if

        nkeep = count(active)
        if (nkeep == 0) then
            call mark_minimizers(0.0_dp, aa_s, bb_s, cc_s, active)
            nkeep = count(active)
        end if

        keep(1:nkeep) = pack(keep, active)
        if (nkeep < size(keep)) keep = keep(1:nkeep)

        deallocate(active, roots, aa_s, bb_s, cc_s)
    end subroutine keep_lower_envelope

    logical pure function quadratic_dominated(a1, b1, c1, a2, b2, c2)
        real(kind=dp), intent(in) :: a1, b1, c1, a2, b2, c2
        real(kind=dp) :: aa, bb, cc, disc

        aa = a1 - a2
        bb = b1 - b2
        cc = c1 - c2
        if (aa < -1.0e-12_dp) then
            quadratic_dominated = .false.
        else if (abs(aa) <= 1.0e-12_dp) then
            if (abs(bb) <= 1.0e-12_dp) then
                quadratic_dominated = (cc >= -1.0e-12_dp)
            else
                quadratic_dominated = .false.
            end if
        else
            disc = bb * bb - 4.0_dp * aa * cc
            quadratic_dominated = (disc <= 1.0e-10_dp)
        end if
    end function quadratic_dominated

    subroutine mark_minimizers(x, a, b, c, active)
        real(kind=dp), intent(in) :: x
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        logical, intent(inout) :: active(:)
        integer :: i
        real(kind=dp) :: minv, v, tol

        minv = a(1) * x * x + b(1) * x + c(1)
        do i = 2, size(a)
            v = a(i) * x * x + b(i) * x + c(i)
            if (v < minv) minv = v
        end do
        tol = 1.0e-8_dp * max(1.0_dp, abs(minv))
        do i = 1, size(a)
            v = a(i) * x * x + b(i) * x + c(i)
            if (v <= minv + tol) active(i) = .true.
        end do
    end subroutine mark_minimizers

    subroutine sort_real(x)
        real(kind=dp), intent(inout) :: x(:)
        integer :: i, j
        real(kind=dp) :: tmp

        do i = 2, size(x)
            tmp = x(i)
            j = i - 1
            do while (j >= 1 .and. x(j) > tmp)
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = tmp
        end do
    end subroutine sort_real

    subroutine accumulate_normal_equations(t, yval, cps, xtx, xty)
        real(kind=dp), intent(in) :: t, yval
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(inout) :: xtx(:, :), xty(:)

        integer :: p, i, j
        real(kind=dp), allocatable :: x(:)

        p = size(cps) + 2
        allocate(x(p))
        x(1) = 1.0_dp
        x(2) = t
        do i = 1, size(cps)
            x(i + 2) = max(0.0_dp, t - real(cps(i), dp))
        end do

        do i = 1, p
            xty(i) = xty(i) + x(i) * yval
            do j = 1, p
                xtx(i, j) = xtx(i, j) + x(i) * x(j)
            end do
        end do

        deallocate(x)
    end subroutine accumulate_normal_equations

    subroutine solve_linear_system(a, b, x)
        real(kind=dp), intent(in) :: a(:, :), b(:)
        real(kind=dp), intent(out) :: x(:)

        integer :: n, i, k, ipiv
        real(kind=dp), allocatable :: aa(:, :), bb(:), rowtmp(:)
        real(kind=dp) :: piv, factor, best

        n = size(b)
        allocate(aa(n, n), bb(n), rowtmp(n))
        aa = a
        bb = b

        do k = 1, n - 1
            ipiv = k
            best = abs(aa(k, k))
            do i = k + 1, n
                if (abs(aa(i, k)) > best) then
                    best = abs(aa(i, k))
                    ipiv = i
                end if
            end do
            if (best <= 1.0e-12_dp) error stop 'solve_linear_system: singular matrix'

            if (ipiv /= k) then
                rowtmp = aa(k, :)
                aa(k, :) = aa(ipiv, :)
                aa(ipiv, :) = rowtmp
                piv = bb(k)
                bb(k) = bb(ipiv)
                bb(ipiv) = piv
            end if

            do i = k + 1, n
                factor = aa(i, k) / aa(k, k)
                aa(i, k:n) = aa(i, k:n) - factor * aa(k, k:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do
        if (abs(aa(n, n)) <= 1.0e-12_dp) error stop 'solve_linear_system: singular matrix'

        x(n) = bb(n) / aa(n, n)
        do i = n - 1, 1, -1
            x(i) = (bb(i) - sum(aa(i, i+1:n) * x(i+1:n))) / aa(i, i)
        end do

        deallocate(aa, bb, rowtmp)
    end subroutine solve_linear_system

end module changepoint_mod
