module io_utils_mod
    use kind_mod, only: dp
    use basic_stats_mod, only: mean, cov_mat, biased_cov_sd
    use util_mod, only: sort_int, cumul_sum
    use pca_jacobi_mod, only: principal_components_cov
    use dataframe_index_date_mod, only: DataFrame_index_date, nrow
    implicit none
    public :: print_true_params, print_model_selection, print_estimated_parameters, &
              print_estimated_parameters_dates, print_corrmat_segment, print_corrmat_diff, &
              corrmat_model_range, print_corrmat_model, keep_obs, &
              print_univar_segments, print_covmat_model, print_pca_loadings, &
              print_lower_corr_sd, print_covmat_segment

contains

    !> Prints the true simulation parameters and empirical correlations for segments.
    subroutine print_true_params(true_cps, corr_true, x, y)
        integer, intent(in) :: true_cps(:)        ! Changepoint locations
        real(kind=dp), intent(in) :: corr_true(:) ! True correlation per segment
        real(kind=dp), intent(in) :: x(:), y(:)   ! Input series (n)
        integer :: n, k, nseg_true, seg_start_true, seg_end_true
        real(kind=dp) :: r_est

        if (size(x) /= size(y)) &
            error stop "print_true_params: size(x) /= size(y)"
        if (size(corr_true) /= size(true_cps) + 1) &
            error stop "print_true_params: size(corr_true) /= size(true_cps) + 1"

        n = size(x)
        nseg_true = size(corr_true)
        print "(a)", "TRUE PARAMETERS"
        print "(a10,a10,a15,a15)", "Start", "End", "Corr_true", "Corr_sim"
        seg_start_true = 1
        do k = 1, nseg_true
            if (k < nseg_true) then
                seg_end_true = true_cps(k)
            else
                seg_end_true = n
            end if
            r_est = sum((x(seg_start_true:seg_end_true)-mean(x(seg_start_true:seg_end_true)))*(y(seg_start_true:seg_end_true)-mean(y(seg_start_true:seg_end_true)))) / &
                    (sqrt(sum((x(seg_start_true:seg_end_true)-mean(x(seg_start_true:seg_end_true)))**2) * sum((y(seg_start_true:seg_end_true)-mean(y(seg_start_true:seg_end_true)))**2)))
            print "(i10,i10,f15.4,f15.4)", seg_start_true, seg_end_true, corr_true(k), r_est
            seg_start_true = seg_end_true + 1
        end do
    end subroutine print_true_params

    !> Prints AIC/BIC statistics for different numbers of segments.
    subroutine print_model_selection(dp_table, parent, best_aic_cp, best_bic_cp, params_per_seg, &
        print_each)
        real(kind=dp), intent(in) :: dp_table(:, :) ! Dynamic programming table of costs (n, max_m)
        integer, intent(in) :: parent(:, :)        ! DP parent pointers (n, max_m)
        integer, intent(out), optional :: best_aic_cp  ! # changepoints chosen by AIC
        integer, intent(out), optional :: best_bic_cp  ! # changepoints chosen by BIC
        integer, intent(in), optional :: params_per_seg ! parameters per segment (default 2)
        logical, intent(in), optional :: print_each
        logical                       :: print_each_
        ! params_per_seg=2: correlation model  → k = 2*m-1 (1 rho/seg + 1 location/break)
        ! params_per_seg=3: mean-shift model   → k = 3*m-1 (mu+sigma^2/seg + 1 location/break)
        integer :: n, max_m, m, cp, k, k_params, pps
        real(kind=dp) :: ll, aic(size(dp_table, 2)), bic(size(dp_table, 2))
        integer :: cps(size(dp_table, 2))
        integer :: best_m_aic, best_m_bic
        real(kind=dp) :: min_aic, min_bic
        print_each_ = .true.
        if (present(print_each)) print_each_ = print_each
        pps = 2
        if (present(params_per_seg)) pps = params_per_seg

        if (any(shape(dp_table) /= shape(parent))) &
            error stop "print_model_selection: dp_table and parent shapes are incompatible"

        n = size(dp_table, 1)
        max_m = size(dp_table, 2)

        if (print_each_) print "(/,a)", "M      LL         AIC         BIC    ChangePoints"
        do m = 1, max_m
            if (dp_table(n, m) >= 1.0e19_dp) then
                aic(m) = 1.0e20_dp
                bic(m) = 1.0e20_dp
                cycle
            end if
            ll = -dp_table(n, m)
            k_params = pps * m - 1
            aic(m) = -2.0_dp*ll + 2.0_dp * k_params
            bic(m) = -2.0_dp*ll + real(k_params, dp) * log(real(n, dp))

            if (print_each_) write(*, "(i2, 3f12.2, 4x)", advance='no') &
               m, ll, aic(m), bic(m)
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))
            if (print_each_) then
               do k = 1, m-1
                  write(*, "(i0, ' ')", advance='no') cps(k)
               end do
               print *
            end if
        end do
        
        min_aic = 1.0e20_dp
        min_bic = 1.0e20_dp
        best_m_aic = 1
        best_m_bic = 1
        do m = 1, max_m
            if (aic(m) < min_aic) then
                min_aic = aic(m)
                best_m_aic = m
            end if
            if (bic(m) < min_bic) then
                min_bic = bic(m)
                best_m_bic = m
            end if
        end do
        print "(/,a,2(1x,i0),/)", "Changepoints chosen by AIC, BIC:", best_m_aic - 1, best_m_bic - 1
        if (present(best_aic_cp)) best_aic_cp = best_m_aic - 1
        if (present(best_bic_cp)) best_bic_cp = best_m_bic - 1
    end subroutine print_model_selection

    !> Prints the estimated parameters (breakpoints and correlations) for each model using indices.
    subroutine print_estimated_parameters(max_m, parent, x, y)
        integer, intent(in) :: max_m              ! Max changepoints
        integer, intent(in) :: parent(:, :)       ! DP parent pointers
        real(kind=dp), intent(in) :: x(:), y(:)   ! Input series
        integer :: n, m, cp, k, seg_start, seg_end
        integer :: cps(max_m)
        real(kind=dp) :: r_est

        if (size(x) /= size(y)) &
            error stop "print_estimated_parameters: size(x) /= size(y)"

        n = size(parent, 1)
        do m = 1, max_m
            ! Reconstruct and sort
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))

            print "(a,i2)", "estimated parameters for m =", m
            print "(a10,a10,a10)", "Start", "End", "Corr"
            seg_start = 1
            do k = 1, m
                if (k == m) then
                    seg_end = n
                else
                    seg_end = cps(k)
                end if
                r_est = sum((x(seg_start:seg_end)-mean(x(seg_start:seg_end)))*(y(seg_start:seg_end)-mean(y(seg_start:seg_end)))) / &
                        (sqrt(sum((x(seg_start:seg_end)-mean(x(seg_start:seg_end)))**2) * sum((y(seg_start:seg_end)-mean(y(seg_start:seg_end)))**2)))
                print "(i10,i10,f10.4)", seg_start, seg_end, r_est
                seg_start = seg_end + 1
            end do
            print *
        end do
    end subroutine print_estimated_parameters

    !> Prints the estimated parameters (breakpoints and correlations) for each model using date labels.
    subroutine print_estimated_parameters_dates(max_m, parent, x, y, dates)
        integer, intent(in) :: max_m              ! Max changepoints
        integer, intent(in) :: parent(:, :)       ! DP parent pointers
        real(kind=dp), intent(in) :: x(:), y(:)   ! Input series
        character(len=*), intent(in) :: dates(:)  ! Date labels
        integer :: n, m, cp, k, seg_start, seg_end
        integer :: cps(max_m)
        real(kind=dp) :: r_est, mx, my, sdx, sdy, covar
        real(kind=dp), allocatable :: x_seg(:), y_seg(:)

        if (size(x) /= size(y)) &
            error stop "print_estimated_parameters_dates: size(x) /= size(y)"

        n = size(parent, 1)
        allocate(x_seg(n), y_seg(n))
        do m = 1, max_m
            ! Reconstruct and sort
            cp = n
            cps(m) = n
            do k = m, 2, -1
                cps(k-1) = parent(cp, k)
                cp = cps(k-1)
            end do
            if (m > 1) call sort_int(cps(1:m-1))

            print "(a,i2)", "estimated parameters for m =", m
            print "(a12,a12,a10,a10,a10,a10,a10,a10,a10)", "Start", "End", "#obs", "Corr", "Covar", "sd_x", "sd_y", "mean_x", "mean_y"
            seg_start = 1
            do k = 1, m
                if (k == m) then
                    seg_end = n
                else
                    seg_end = cps(k)
                end if
                
                ! Extract segment
                x_seg(1:seg_end-seg_start+1) = x(seg_start:seg_end)
                y_seg(1:seg_end-seg_start+1) = y(seg_start:seg_end)
                mx = mean(x_seg(1:seg_end-seg_start+1))
                my = mean(y_seg(1:seg_end-seg_start+1))
                sdx = sqrt(sum((x_seg(1:seg_end-seg_start+1)-mx)**2) / max(1, (seg_end-seg_start)))
                sdy = sqrt(sum((y_seg(1:seg_end-seg_start+1)-my)**2) / max(1, (seg_end-seg_start)))
                covar = sum((x_seg(1:seg_end-seg_start+1)-mx)*(y_seg(1:seg_end-seg_start+1)-my)) / max(1, (seg_end-seg_start))
                r_est = covar / (sdx * sdy)
                
                ! r_est, covar, sdx, sdy, mx, my
                print "(a12,a12,i10,f10.4,5f10.4)", dates(seg_start), dates(seg_end), seg_end - seg_start + 1, r_est, covar, sdx, sdy, mx, my

                seg_start = seg_end + 1
            end do
            print *
        end do
        deallocate(x_seg, y_seg)
    end subroutine print_estimated_parameters_dates
    !> Prints correlation matrix (lower triangle), annualised std devs, and
    !! annualised arithmetic returns for one segment of a corrmat changepoint model.
    subroutine print_corrmat_segment(k, i0, i1, R, col_names, dates, p, scale_ret, in_sd, in_mu, do_pca_cov, do_pca_corr)
        integer,          intent(in)           :: k, i0, i1, p
        real(kind=dp),    intent(in)           :: R(:,:)
        character(len=*), intent(in)           :: col_names(:), dates(:)
        real(kind=dp),    intent(in)           :: scale_ret
        real(kind=dp),    intent(in), optional :: in_sd(:), in_mu(:)
        logical,          intent(in), optional :: do_pca_cov, do_pca_corr
        integer       :: a, b, m
        real(kind=dp) :: mu(p), sd(p), cov_ab, r_ab
        real(kind=dp) :: S(p, p), C(p, p), sd_s(p)
        real(kind=dp), parameter :: ann = 15.87401_dp  ! sqrt(252), daily → annual

        m = i1 - i0 + 1
        print "(/,'Segment ',i0,': ',a,' to ',a,' (',i0,' obs)')", &
            k, trim(dates(i0)), trim(dates(i1)), m

        if (present(in_sd) .and. present(in_mu)) then
            sd = in_sd
            mu = in_mu
        else
            do a = 1, p
                mu(a) = sum(R(i0:i1, a)) / m
                sd(a) = sqrt(max(sum((R(i0:i1,a) - mu(a))**2) / m, 0.0_dp))
            end do
        end if

        print "(a)", "  Correlation:"
        write (*, "(8x)", advance="no")
        do a = 1, p
            write (*, "(a8)", advance="no") trim(col_names(a))
        end do
        print *
        do a = 1, p
            write (*, "(4x,a4)", advance="no") trim(col_names(a))
            do b = 1, a
                if (sd(a) > 0.0_dp .and. sd(b) > 0.0_dp) then
                    cov_ab = sum((R(i0:i1,a) - mu(a)) * (R(i0:i1,b) - mu(b))) / m
                    r_ab   = cov_ab / (sd(a) * sd(b))
                else
                    r_ab = 0.0_dp
                end if
                write (*, "(f8.3)", advance="no") r_ab
            end do
            print *
        end do

        write (*, "(4x,a4,*(f8.3))") "*SD*",  sd * ann / scale_ret
        write (*, "(3x,a5,*(f8.3))") "*RET*", mu * 252 / scale_ret

        if ((present(do_pca_cov) .and. do_pca_cov) .or. &
            (present(do_pca_corr) .and. do_pca_corr)) then
            call biased_cov_sd(R(i0:i1, :), S, sd_s)
            if (present(do_pca_cov) .and. do_pca_cov) &
                call print_pca_loadings(S, col_names, "PCA of covariance matrix")
            if (present(do_pca_corr) .and. do_pca_corr) then
                do a = 1, p
                    do b = 1, p
                        if (sd_s(a) > 0.0_dp .and. sd_s(b) > 0.0_dp) then
                            C(a,b) = S(a,b) / (sd_s(a) * sd_s(b))
                        else
                            C(a,b) = merge(1.0_dp, 0.0_dp, a == b)
                        end if
                    end do
                end do
                call print_pca_loadings(C, col_names, "PCA of correlation matrix")
            end if
        end if
    end subroutine print_corrmat_segment

    !> For each pair of assets, tests whether the correlation changed significantly
    !! between two segments using the Fisher z-test with Bonferroni correction.
    !! Only pairs that are significant after correction are printed.
    !!
    !! H0: rho1 = rho2.  Test statistic: (atanh(r1)-atanh(r2)) / sqrt(1/(n1-3)+1/(n2-3))
    !! Two-sided p-value via erfc.  Bonferroni threshold: alpha / (p*(p-1)/2).
    subroutine print_corrmat_diff(i0a, i1a, i0b, i1b, R, col_names, alpha, sd1, mu1, sd2, mu2, scale_ret)
        integer,          intent(in) :: i0a, i1a, i0b, i1b
        real(kind=dp),    intent(in) :: R(:,:)
        character(len=*), intent(in) :: col_names(:)
        real(kind=dp),    intent(in) :: alpha, scale_ret
        real(kind=dp),    intent(in) :: sd1(:), mu1(:), sd2(:), mu2(:)

        integer       :: p, ia, ib, na, nb, npairs, nsig, k
        real(kind=dp) :: r1, r2, z1, z2, se, stat, p_val, thresh
        real(kind=dp), parameter :: r_clamp = 1.0_dp - 1.0e-10_dp
        real(kind=dp), parameter :: ann = 15.87401_dp  ! sqrt(252), daily → annual
        integer,       allocatable :: sig_ia(:), sig_ib(:)
        real(kind=dp), allocatable :: sig_r1(:), sig_r2(:), sig_stat(:), sig_pval(:)
        character(len=11) :: pair_str

        p      = size(col_names)
        na     = i1a - i0a + 1
        nb     = i1b - i0b + 1
        npairs = p * (p - 1) / 2
        thresh = alpha / npairs

        allocate(sig_ia(npairs), sig_ib(npairs), &
                 sig_r1(npairs), sig_r2(npairs), sig_stat(npairs), sig_pval(npairs))

        nsig = 0
        do ia = 1, p
            do ib = ia + 1, p
                if (sd1(ia) > 0.0_dp .and. sd1(ib) > 0.0_dp) then
                    r1 = sum((R(i0a:i1a,ia)-mu1(ia)) * (R(i0a:i1a,ib)-mu1(ib))) &
                         / (na * sd1(ia) * sd1(ib))
                else
                    r1 = 0.0_dp
                end if
                if (sd2(ia) > 0.0_dp .and. sd2(ib) > 0.0_dp) then
                    r2 = sum((R(i0b:i1b,ia)-mu2(ia)) * (R(i0b:i1b,ib)-mu2(ib))) &
                         / (nb * sd2(ia) * sd2(ib))
                else
                    r2 = 0.0_dp
                end if
                r1    = max(-r_clamp, min(r_clamp, r1))
                r2    = max(-r_clamp, min(r_clamp, r2))
                z1    = atanh(r1)
                z2    = atanh(r2)
                se    = sqrt(1.0_dp/(na - 3) + 1.0_dp/(nb - 3))
                stat  = (z1 - z2) / se
                p_val = erfc(abs(stat) / sqrt(2.0_dp))
                if (p_val < thresh) then
                    nsig = nsig + 1
                    sig_ia(nsig)   = ia
                    sig_ib(nsig)   = ib
                    sig_r1(nsig)   = r1
                    sig_r2(nsig)   = r2
                    sig_stat(nsig) = stat
                    sig_pval(nsig) = p_val
                end if
            end do
        end do

        if (nsig > 0) then
            print "('  Sig. corr. changes (Bonferroni adj. alpha=',f6.4,'): ',i0,'/',i0,' = ',f5.3)", &
                thresh, nsig, npairs, real(nsig, dp) / npairs
            print "(4x, a11, 4a8, a9, 8a7)", "Pair       ", "r1", "r2", "r1-r2", "z", "p", &
                "sd_x1", "sd_x2", "sd_y1", "sd_y2", "ret_x1", "ret_x2", "ret_y1", "ret_y2"
            do k = 1, nsig
                pair_str = trim(col_names(sig_ia(k))) // "-" // trim(col_names(sig_ib(k)))
                print "(4x, a11, 4f8.3, f9.6, 8f7.3)", pair_str, &
                    sig_r1(k), sig_r2(k), sig_r1(k) - sig_r2(k), sig_stat(k), sig_pval(k), &
                    sd1(sig_ia(k)) * ann / scale_ret, sd2(sig_ia(k)) * ann / scale_ret, &
                    sd1(sig_ib(k)) * ann / scale_ret, sd2(sig_ib(k)) * ann / scale_ret, &
                    mu1(sig_ia(k)) * 252 / scale_ret, mu2(sig_ia(k)) * 252 / scale_ret, &
                    mu1(sig_ib(k)) * 252 / scale_ret, mu2(sig_ib(k)) * 252 / scale_ret
            end do
        end if

        deallocate(sig_ia, sig_ib, sig_r1, sig_r2, sig_stat, sig_pval)
    end subroutine print_corrmat_diff

    subroutine print_corrmat_model(best_bic, seg_ends, R, col_names, ret_dates, &
                                   scale_ret, print_diffs, alpha_diff, do_pca_cov, do_pca_corr)
        !> Print the correlation structure for one changepoint model.
        !! Pre-computes segment means and std devs once and passes them to
        !! print_corrmat_segment and print_corrmat_diff to avoid recomputation.
        integer,          intent(in)           :: best_bic, seg_ends(:)
        real(kind=dp),    intent(in)           :: R(:,:)
        character(len=*), intent(in)           :: col_names(:), ret_dates(:)
        real(kind=dp),    intent(in)           :: scale_ret, alpha_diff
        logical,          intent(in)           :: print_diffs
        logical,          intent(in), optional :: do_pca_cov, do_pca_corr
        integer :: m_segs, n_col, k, seg_start, m_k, a
        real(kind=dp), allocatable :: seg_mu(:,:), seg_sd(:,:)

        m_segs = size(seg_ends)
        n_col  = size(col_names)

        if (m_segs - 1 == best_bic) then
            print "(/,'Correlation structure: ',i0,' changepoint',a,' (',i0,' segment',a,') [BIC]')", &
                m_segs-1, merge("s"," ", m_segs-1 /= 1), m_segs, merge("s"," ", m_segs /= 1)
        else
            print "(/,'Correlation structure: ',i0,' changepoint',a,' (',i0,' segment',a,')')", &
                m_segs-1, merge("s"," ", m_segs-1 /= 1), m_segs, merge("s"," ", m_segs /= 1)
        end if

        allocate(seg_mu(n_col, m_segs), seg_sd(n_col, m_segs))
        seg_start = 1
        do k = 1, m_segs
            m_k = seg_ends(k) - seg_start + 1
            do a = 1, n_col
                seg_mu(a, k) = sum(R(seg_start:seg_ends(k), a)) / m_k
                seg_sd(a, k) = sqrt(max(sum((R(seg_start:seg_ends(k), a) - seg_mu(a,k))**2) / m_k, 0.0_dp))
            end do
            seg_start = seg_ends(k) + 1
        end do

        seg_start = 1
        do k = 1, m_segs
            call print_corrmat_segment(k, seg_start, seg_ends(k), R, col_names, ret_dates, n_col, &
                                       scale_ret, in_sd=seg_sd(:,k), in_mu=seg_mu(:,k), &
                                       do_pca_cov=do_pca_cov, do_pca_corr=do_pca_corr)
            if (print_diffs .and. k < m_segs) &
                call print_corrmat_diff(seg_start, seg_ends(k), seg_ends(k)+1, seg_ends(k+1), &
                                        R, col_names, alpha_diff, &
                                        seg_sd(:,k), seg_mu(:,k), seg_sd(:,k+1), seg_mu(:,k+1), scale_ret)
            seg_start = seg_ends(k) + 1
        end do
        deallocate(seg_mu, seg_sd)
    end subroutine print_corrmat_model

    pure elemental subroutine corrmat_model_range(print_segs, best_bic, max_m, m_lo, m_hi)
        !> Maps the print_segs control parameter to a range [m_lo, m_hi] of segment
        !! counts to display: 0 = BIC-chosen only, 1 = 0 through BIC, 2 = all studied.
        integer, intent(in)  :: print_segs, best_bic, max_m
        integer, intent(out) :: m_lo, m_hi
        select case (print_segs)
        case (1)
            m_lo = 1
            m_hi = best_bic + 1
        case (2)
            m_lo = 1
            m_hi = max_m
        case default  ! 0: BIC only
            m_lo = best_bic + 1
            m_hi = best_bic + 1
        end select
    end subroutine corrmat_model_range

    subroutine keep_obs(df, max_days, latest, verbose)
        !> Keep at most max_days rows of df; print the range used if verbose is present and true.
        type(DataFrame_index_date), intent(in out) :: df
        integer, intent(in) :: max_days
        logical, intent(in) :: latest
        logical, intent(in), optional :: verbose
        df = df%keep_rows(max_days + 1, latest=latest)
        if (present(verbose)) then
            if (verbose) print "('using ',a,' ',i0,' values (',a,' to ',a,')')", &
                merge("latest  ", "earliest", latest), nrow(df) - 1, &
                trim(df%index(2)%to_str()), trim(df%index(nrow(df))%to_str())
        end if
    end subroutine keep_obs

    subroutine print_univar_segments(best_bic, seg_ends, z, series_name, ret_dates)
        !> Print segment means and standard deviations for a univariate series z.
        integer, intent(in)          :: best_bic, seg_ends(:)
        real(kind=dp), intent(in)    :: z(:)
        character(len=*), intent(in) :: series_name, ret_dates(:)
        integer :: ms, k, i0, i1, nseg
        real(kind=dp) :: zmean, zsd
        ms = size(seg_ends)
        if (ms == best_bic + 1) then
            print "(/,'BIC-selected model (',i0,' changepoint(s)) -- ',a)", best_bic, series_name
        else
            print "(/,'model (',i0,' changepoint(s)) -- ',a)", ms - 1, series_name
        end if
        print "(a12,a12,a8,a12,a12)", "start", "end", "n", "mean", "sd"
        i0 = 1
        do k = 1, ms
            i1 = seg_ends(k)
            nseg = i1 - i0 + 1
            zmean = sum(z(i0:i1)) / nseg
            if (nseg > 1) then
                zsd = sqrt(sum((z(i0:i1) - zmean)**2) / (nseg - 1))
            else
                zsd = 0.0_dp
            end if
            print "(a12,a12,i8,2f12.4)", ret_dates(i0), ret_dates(i1), nseg, zmean, zsd
            i0 = i1 + 1
        end do
    end subroutine print_univar_segments

    subroutine print_covmat_model(best_bic, seg_ends, R, col_names, ret_dates)
        !> Print annualized covariance matrix (x252) for each segment.
        integer, intent(in)          :: best_bic, seg_ends(:)
        real(kind=dp), intent(in)    :: R(:,:)
        character(len=*), intent(in) :: col_names(:), ret_dates(:)
        real(kind=dp), parameter     :: days_per_year = 252.0_dp
        integer :: ms, k, jc, i0, i1
        real(kind=dp), allocatable :: C(:,:)
        ms = size(seg_ends)
        if (ms == best_bic + 1) then
            print "(/,'BIC-selected model (',i0,' changepoint(s))')", best_bic
        else
            print "(/,'model (',i0,' changepoint(s))')", ms - 1
        end if
        i0 = 1
        do k = 1, ms
            i1 = seg_ends(k)
            print "('segment ',i0,': ',a,' to ',a,' (',i0,' obs)')", &
                k, ret_dates(i0), ret_dates(i1), i1 - i0 + 1
            C = cov_mat(R(i0:i1, :)) * days_per_year
            print "(5x,*(a8,:,1x))", (trim(col_names(jc)), jc=1,size(col_names))
            do jc = 1, size(col_names)
                print "(a8,*(1x,f8.4))", trim(col_names(jc)), C(jc,:)
            end do
            i0 = i1 + 1
        end do
    end subroutine print_covmat_model

    subroutine print_covmat_segment(k, i0, i1, R, col_names, dates, scale_ret, do_pca)
        !> Print one segment of a covariance changepoint model: header, correlation
        !! lower triangle, annualized std devs, and optionally PCA loadings.
        integer,          intent(in) :: k, i0, i1
        real(kind=dp),    intent(in) :: R(:,:)
        character(len=*), intent(in) :: col_names(:), dates(:)
        real(kind=dp),    intent(in) :: scale_ret
        logical,          intent(in) :: do_pca
        real(kind=dp) :: S(size(col_names), size(col_names)), sd(size(col_names))
        print "(/,'Segment ',i0,': ',a,' to ',a,' (',i0,' obs)')", &
            k, trim(dates(i0)), trim(dates(i1)), i1 - i0 + 1
        call biased_cov_sd(R(i0:i1, :), S, sd)
        call print_lower_corr_sd(S, sd, col_names, scale_ret)
        if (do_pca) call print_pca_loadings(S, col_names)
    end subroutine print_covmat_segment

    subroutine print_lower_corr_sd(S, sd, col_names, scale_ret)
        !> Print correlation matrix (lower triangle) derived from covariance matrix S
        !! and standard deviations sd, followed by a row of annualized std devs.
        real(kind=dp),    intent(in) :: S(:,:), sd(:)
        character(len=*), intent(in) :: col_names(:)
        real(kind=dp),    intent(in) :: scale_ret
        real(kind=dp), parameter :: ann = 15.87401_dp  ! sqrt(252), daily -> annual
        integer :: a, b, p
        real(kind=dp) :: r_ab
        p = size(col_names)
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
        write(*, "(4x,a4)", advance='no') '*SD*'
        do a = 1, p
            write(*, "(f8.3)", advance='no') sd(a) * ann / scale_ret
        end do
        print *
    end subroutine print_lower_corr_sd

    subroutine print_pca_loadings(S, col_names, title)
        !> Print principal component loadings, variance explained, and cumulative
        !! variance explained for matrix S with given column labels and optional title.
        real(kind=dp),    intent(in)           :: S(:,:)
        character(len=*), intent(in)           :: col_names(:)
        character(len=*), intent(in), optional :: title
        real(kind=dp), allocatable :: evals(:), evecs(:,:), var_explained(:)
        integer :: a, j, p
        p = size(col_names)
        call principal_components_cov(S, evals, evecs, var_explained)
        if (present(title)) then
            print "(/,a)", "  " // title // ":"
        else
            print "(/,a)", "  Principal component loadings:"
        end if
        write(*, "(a12,*(1x,f12.3))") "  var_exp", var_explained
        write(*, "(a12,*(1x,f12.3))") "  cumul",   cumul_sum(var_explained)
        write(*, "(a12,*(1x,i12))")   "  PC",       (j, j=1,p)
        do a = 1, p
            write(*, "(a12,*(1x,f12.3))") "  " // trim(col_names(a)), evecs(a,:)
        end do
    end subroutine print_pca_loadings

end module io_utils_mod
