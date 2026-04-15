program xreturns_pca
! read prices, compute returns and principal components of covariance matrix
use kind_mod, only: dp, long_int
use dataframe_index_date_mod, only: dataframe_index_date, ncol, print_summary, operator(*)
use basic_stats_mod, only: cov_mat
use pca_jacobi_mod, only: principal_components_cov
use util_mod, only: print_wall_time, cumul_sum, print_square_matrix
implicit none

type(dataframe_index_date) :: df_px, df_ret
real(kind=dp), allocatable :: xcov(:,:), evals(:), evecs(:,:), var_explained(:)
integer(kind=long_int) :: t_start
integer :: i, j, p
real(kind=dp), parameter :: scale_ret = 100.0_dp
character(len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv"

call system_clock(t_start)
call df_px%read_csv(prices_file)
print "(a)", "prices file: " // trim(prices_file)

df_ret = scale_ret * df_px%pct_change()
print "(a, f0.1)", "returns scaled by: ", scale_ret
call print_summary(df_ret)

xcov = cov_mat(df_ret%values(2:, :))
call principal_components_cov(xcov, evals, evecs, var_explained)

p = size(evals)
call print_square_matrix(xcov, df_ret%columns, "covariance matrix")

print "(/,a)", "principal component loadings"
print "(a12,*(1x,f12.3))", "var_exp", var_explained ! (var_explained(j), j=1,p)
print "(a12,*(1x,f12.3))", "cumul", cumul_sum(var_explained) ! (var_explained(j), j=1,p)
print "(a12,*(1x,i12))", "PC", (j, j=1,p)
do i = 1, p
   write (*,"(a12,*(1x,f12.3))") trim(df_ret%columns(i)), evecs(i,:)
end do

call print_wall_time(t_start)
end program xreturns_pca

