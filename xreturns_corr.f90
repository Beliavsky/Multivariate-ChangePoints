program xreturns_corr
use kind_mod, only: dp
use dataframe_index_date_mod, only: DataFrame_index_date, ncol, print_summary, &
   operator(*)
use basic_stats_mod, only: col_stats_ignore_nan, print_corr_mat, print_acf, rms
use util_mod, only: print_time_elapsed, print_wall_time
implicit none

type(DataFrame_index_date) :: df_px, df_ret
integer :: j, n
integer(kind=8) :: t_start
real(kind=dp) :: mean_ret, sd_ret, min_ret, max_ret
real(kind=dp), parameter :: obs_year = 252.0_dp, scale_ret = 100.0_dp
logical, parameter       :: use_rms = .true.         ! if .true., use RMS (zero mean) instead of SD for volatilities
logical, parameter       :: print_corr_ret = .true.
character (len=*), parameter :: prices_file = "spy_efa_eem_tlt.csv"
real(kind=dp) :: xvol
call system_clock(t_start)
call df_px%read_csv(prices_file) ! read prices from csv file
print "(a)", "prices file: " // trim(prices_file)
df_ret = scale_ret * df_px%pct_change() ! compute simple returns
print "(a, f0.1)", "returns scaled by: ", scale_ret
print "(a, l1)", "use_rms: ", use_rms
call print_summary(df_ret) ! print summary of returns dataframe

print "(/,a)", "return statistics by column"
print "(a12,1x,a8,*(1x,a8))", "column", "n","min","max","mean_ann","sd_ann"

do j = 1, ncol(df_ret)
   call col_stats_ignore_nan(df_ret%values(:,j), n, mean_ret, sd_ret, min_ret, max_ret) ! compute stats ignoring nan values
   if (use_rms) then 
      xvol = sqrt(obs_year) * rms(df_ret%values(2:,j))
   else
      xvol = sqrt(obs_year) * sd_ret
   end if
   write (*,"(a12,1x,i8,4(1x,f8.2))") trim(df_ret%columns(j)), n, &
      min_ret, max_ret, obs_year*mean_ret, xvol
end do

if (print_corr_ret) then
   print "(/,a)", "correlation matrix of returns"
   call print_corr_mat(df_ret%values(2:, :), df_ret%columns)
end if
call print_wall_time(t_start)
end program xreturns_corr
