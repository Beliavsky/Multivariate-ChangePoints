module pca_jacobi_mod
use kind_mod, only: dp
implicit none
private
public :: principal_components_cov, jacobi_eigen_sym

contains

subroutine principal_components_cov(covmat, evals, evecs, var_explained, tol, max_sweeps)
! compute principal components of a symmetric covariance matrix
real(kind=dp), intent(in) :: covmat(:,:)
real(kind=dp), allocatable, intent(out) :: evals(:)
real(kind=dp), allocatable, intent(out) :: evecs(:,:)
real(kind=dp), allocatable, intent(out), optional :: var_explained(:)
real(kind=dp), intent(in), optional :: tol
integer, intent(in), optional :: max_sweeps
real(kind=dp) :: total_var
integer :: n

n = size(covmat, 1)
if (size(covmat, 2) /= n) then
   error stop "principal_components_cov: covmat must be square"
end if
if (n <= 0) then
   error stop "principal_components_cov: covmat must have positive size"
end if
if (maxval(abs(covmat - transpose(covmat))) > 100.0_dp*epsilon(1.0_dp)) then
   error stop "principal_components_cov: covmat must be symmetric"
end if

call jacobi_eigen_sym(covmat, evals, evecs, tol, max_sweeps)

if (present(var_explained)) then
   allocate(var_explained(n))
   total_var = sum(evals)
   if (total_var > 0.0_dp) then
      var_explained = evals / total_var
   else
      var_explained = 0.0_dp
   end if
end if
end subroutine principal_components_cov

subroutine jacobi_eigen_sym(a_in, evals, evecs, tol, max_sweeps)
! compute eigenpairs of a real symmetric matrix by Jacobi rotations
real(kind=dp), intent(in) :: a_in(:,:)
real(kind=dp), allocatable, intent(out) :: evals(:)
real(kind=dp), allocatable, intent(out) :: evecs(:,:)
real(kind=dp), intent(in), optional :: tol
integer, intent(in), optional :: max_sweeps
real(kind=dp), allocatable :: a(:,:), v(:,:), tmp(:)
real(kind=dp) :: tol_, app, aqq, apq, c, s, tau, t, offmax
real(kind=dp) :: aip, aiq, vip, viq
integer :: n, i, j, p, q, sweep, max_sweeps_
logical :: converged

n = size(a_in, 1)
if (size(a_in, 2) /= n) then
   error stop "jacobi_eigen_sym: matrix must be square"
end if
if (maxval(abs(a_in - transpose(a_in))) > 100.0_dp*epsilon(1.0_dp)) then
   error stop "jacobi_eigen_sym: matrix must be symmetric"
end if

allocate(a(n,n), v(n,n), evals(n), evecs(n,n), tmp(n))
a = a_in
v = 0.0_dp
do i = 1, n
   v(i,i) = 1.0_dp
end do

tol_ = sqrt(epsilon(1.0_dp))
if (present(tol)) tol_ = tol
max_sweeps_ = max(20*n*n, 50)
if (present(max_sweeps)) max_sweeps_ = max_sweeps

converged = .false.
do sweep = 1, max_sweeps_
   offmax = 0.0_dp
   p = 1
   q = 1
   do j = 2, n
      do i = 1, j - 1
         if (abs(a(i,j)) > offmax) then
            offmax = abs(a(i,j))
            p = i
            q = j
         end if
      end do
   end do

   if (offmax <= tol_) then
      converged = .true.
      exit
   end if

   app = a(p,p)
   aqq = a(q,q)
   apq = a(p,q)

   if (apq /= 0.0_dp) then
      tau = (aqq - app) / (2.0_dp*apq)
      if (tau >= 0.0_dp) then
         t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
         t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t)
      s = t*c
   else
      c = 1.0_dp
      s = 0.0_dp
   end if

   do i = 1, n
      if (i /= p .and. i /= q) then
         aip = a(i,p)
         aiq = a(i,q)
         a(i,p) = c*aip - s*aiq
         a(p,i) = a(i,p)
         a(i,q) = s*aip + c*aiq
         a(q,i) = a(i,q)
      end if
   end do

   a(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
   a(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
   a(p,q) = 0.0_dp
   a(q,p) = 0.0_dp

   do i = 1, n
      vip = v(i,p)
      viq = v(i,q)
      v(i,p) = c*vip - s*viq
      v(i,q) = s*vip + c*viq
   end do
end do

if (.not. converged) then
   error stop "jacobi_eigen_sym: no convergence"
end if

do i = 1, n
   evals(i) = a(i,i)
end do
evecs = v
call sort_eigenpairs_desc(evals, evecs)
end subroutine jacobi_eigen_sym

subroutine sort_eigenpairs_desc(evals, evecs)
! sort eigenvalues descending and permute eigenvectors to match
real(kind=dp), intent(in out) :: evals(:)
real(kind=dp), intent(in out) :: evecs(:,:)
real(kind=dp) :: x
real(kind=dp), allocatable :: vtmp(:)
integer :: i, j, k, n

n = size(evals)
allocate(vtmp(size(evecs,1)))
do i = 1, n - 1
   k = i
   do j = i + 1, n
      if (evals(j) > evals(k)) k = j
   end do
   if (k /= i) then
      x = evals(i)
      evals(i) = evals(k)
      evals(k) = x
      vtmp = evecs(:,i)
      evecs(:,i) = evecs(:,k)
      evecs(:,k) = vtmp
   end if
end do
end subroutine sort_eigenpairs_desc

end module pca_jacobi_mod
