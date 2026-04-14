module util_mod
use iso_fortran_env, only: output_unit
use kind_mod, only: dp, long_int
implicit none
private
public :: default, assert_equal, write_merge, split_string, display, &
   print_time_elapsed, print_wall_time, read_words_line, str, print_table, exe_name, &
   join, seq, cbind, sort_int, set_segment_values
interface default
   module procedure default_int, default_real, default_logical, &
      default_character
end interface default
interface seq
   module procedure seq_stride, seq_unit_stride
end interface seq
interface cbind
   module procedure cbind_vec_vec, cbind_mat_vec, cbind_mat_mat
end interface cbind
interface display
   module procedure display_matrix, display_vector
end interface display
contains

    !> Creates a time series by mapping segments to values defined by changepoints.
    function set_segment_values(n, changepoints, values) result(series)
        integer, intent(in) :: n                   ! Length of the time series
        integer, intent(in) :: changepoints(:)     ! Changepoint indices
        real(kind=dp), intent(in) :: values(:)     ! Values per segment
        real(kind=dp) :: series(n)                 ! Resulting time series
        integer :: i, k, j, n_segments

        n_segments = size(values)
        do i = 1, n
            k = n_segments
            do j = 1, n_segments - 1
                if (i <= changepoints(j)) then
                    k = j
                    exit
                end if
            end do
            series(i) = values(k)
        end do
    end function set_segment_values

subroutine sort_int(arr)
    ! Sorts an integer array in ascending order using a simple bubble sort
    integer, intent(in out) :: arr(:)
    integer :: n, i, j, temp
    n = size(arr)
    do i = 1, n - 1
        do j = i + 1, n
            if (arr(i) > arr(j)) then
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
            end if
        end do
    end do
end subroutine sort_int

elemental function default_int(x, xopt) result(y)
integer, intent(in) :: x
integer, intent(in), optional :: xopt
integer             :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_int

elemental function default_real(x, xopt) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp), intent(in), optional :: xopt
real(kind=dp)             :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_real

elemental function default_logical(x, xopt) result(y)
logical, intent(in) :: x
logical, intent(in), optional :: xopt
logical             :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_logical

elemental function default_character(x, xopt) result(y)
character (len=*), intent(in) :: x
character (len=*), intent(in), optional :: xopt
character (len=100) :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_character

subroutine assert_equal(k, kreq, msg)
integer, intent(in) :: k, kreq
character (len=*), intent(in) :: msg
if (k /= kreq) then
   print "(a, i0, a, i0)", msg // " = ", k, ", must equal ", kreq
   stop
end if
end subroutine assert_equal

subroutine write_merge(tf, x, y, outu, fmt)
logical, intent(in) :: tf
character (len=*), intent(in) :: x, y
integer, intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt
integer :: outu_
character (len=100) :: fmt_
outu_ = default(output_unit, outu)
if (present(fmt)) then
   fmt_ = fmt
else
   fmt_ = "(a)"
end if
if (tf) then
   write (outu_, fmt_) x
else
   write (outu_, fmt_) y
end if
end subroutine write_merge

subroutine split_string(str, delim, tokens)
character(len=*), intent(in)           :: str
character(len=*), intent(in)           :: delim
character(:), allocatable, intent(out) :: tokens(:)
integer :: start, pos, i, count, n

n = len_trim(str)
if (n == 0) then
   allocate(character(len=0) :: tokens(1))
   tokens(1) = ""
   return
end if

count = 0
start = 1
do
   pos = index(str(start:), delim)
   if (pos == 0) then
      count = count + 1
      exit
   else
      count = count + 1
      start = start + pos
   end if
end do

allocate(character(len=n) :: tokens(count))

start = 1
i = 1
do
   pos = index(str(start:), delim)
   if (pos == 0) then
      tokens(i) = adjustl(str(start:))
      exit
   else
      tokens(i) = adjustl(str(start:start+pos-2))
      start = start + pos
      i = i + 1
   end if
end do
end subroutine split_string

subroutine display_matrix(x, outu, fmt_r, fmt_header, fmt_trailer, &
   title)
real(kind=dp)    , intent(in)           :: x(:,:)
integer          , intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt_r, fmt_header, &
                                           fmt_trailer, title
integer                                 :: i, outu_
character (len=100)                     :: fmt_r_
outu_  = default(output_unit, outu)
fmt_r_ = default("(*(1x,f10.4))", fmt_r)
if (present(fmt_header)) write(outu_, fmt_header)
if (present(title)) write (outu_, "(a)") title
do i=1,size(x,1)
   write(outu_,fmt_r_) x(i,:)
end do
if (present(fmt_trailer)) write(outu_, fmt_trailer)
end subroutine display_matrix

subroutine display_vector(x, outu, fmt_r, fmt_header, title)
real(kind=dp)    , intent(in)           :: x(:)
integer          , intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt_r, fmt_header, title
integer                                 :: i, outu_
character (len=100)                     :: fmt_r_
outu_  = default(output_unit, outu)
fmt_r_ = default("(*(1x,f10.4))", fmt_r)
if (present(fmt_header)) write(outu_, fmt_header)
if (present(title)) write (outu_, "(a)") title
do i=1,size(x)
   write(outu_,fmt_r_) x(i)
end do
end subroutine display_vector

subroutine print_time_elapsed(old_time, outu)
real(kind=dp), intent(in) :: old_time
real(kind=dp)             :: tt
integer      , intent(in), optional :: outu
integer                             :: outu_
character (len=100) :: fmt_time_
outu_ = default(output_unit, outu)
call cpu_time(tt)
fmt_time_= "('cpu time elapsed (s): ', f0.4)"
write (outu_, fmt_time_) tt - old_time
end subroutine print_time_elapsed

subroutine print_wall_time(t_start, outu)
   integer(kind=long_int), intent(in) :: t_start
   integer(kind=long_int) :: t_end, t_rate
   integer      , intent(in), optional :: outu
   integer                             :: outu_
   outu_ = default(output_unit, outu)
   call system_clock(t_end, t_rate)
   write (outu_, "(/,'wall time elapsed (s): ', f10.4)") real(t_end - t_start, kind=dp) / real(t_rate, kind=dp)
end subroutine print_wall_time


subroutine read_words_line(iu,words)
integer          , intent(in)               :: iu
character (len=*), intent(out), allocatable :: words(:)
integer :: ierr, nwords
character (len=10000) :: text
read (iu,"(a)") text
read (text, *) nwords
allocate (words(nwords))
read (text, *, iostat=ierr) nwords, words
if (ierr /= 0) then
   print*,"could not read ", nwords, " words from '" // trim(text) // "'"
   error stop
end if
end subroutine read_words_line

function str(i) result(text)
integer, intent(in) :: i
character (len=20) :: text
write (text,"(i0)") i
end function str

subroutine print_table(x, row_names, col_names, outu, &
   fmt_col_names, fmt_row, fmt_header, fmt_trailer)
real(kind=dp)    , intent(in) :: x(:,:) 
character (len=*), intent(in) :: row_names(:), col_names(:)
integer          , intent(in), optional :: outu 
character (len=*), intent(in), optional :: fmt_col_names, fmt_row, &
   fmt_header, fmt_trailer
integer                       :: i, n1, n2, outu_
character (len=*), parameter  :: msg="in print_table, "
character (len=100) :: fmt_col_names_, fmt_row_
n1 = size(x, 1)
n2 = size(x, 2)
call assert_equal(size(row_names), n1, msg // "size(row_names)")
call assert_equal(size(col_names), n2, msg // "size(col_names)")
fmt_col_names_ = default("(*(a12,:,1x))", fmt_col_names)
fmt_row_ = default("(a12, *(1x,f12.6))", fmt_row)
outu_ = default(output_unit, outu)
if (present(fmt_header)) write (outu_, fmt_header)
write (outu_, fmt_col_names_) "", (trim(col_names(i)), i=1,n2)
do i=1,n1
   write (outu_, fmt_row_) trim(row_names(i)), x(i,:)
end do
if (present(fmt_trailer)) write (outu_, fmt_trailer)
end subroutine print_table

function exe_name() result(xname)
character (len=1000) :: xname
call get_command_argument(0,xname)
xname = trim(xname)
end function exe_name

function join(words,sep) result(str)
character (len=*), intent(in)                                   :: words(:),sep
character (len=(size(words)-1)*len(sep) + sum(len_trim(words))) :: str
integer                                                         :: i,nw
nw  = size(words)
str = ""
if (nw < 1) then
   return
else
   str = words(1)
end if
do i=2,nw
   str = trim(str) // sep // words(i) 
end do
end function join

pure function seq_stride(first, last, stride) result(vec)
integer, intent(in) :: first, last, stride
integer, allocatable :: vec(:)
integer :: i, n, idiff
idiff = last - first
n = max(0, 1 + idiff/stride)
allocate (vec(n))
do i=1, n
   vec(i) = first + (i - 1) * stride
end do
end function seq_stride

pure function seq_unit_stride(first, last) result(vec)
integer, intent(in) :: first, last
integer, allocatable :: vec(:)
integer :: i, n
n = max(0, last - first + 1)
allocate (vec(n))
do i=1, n
   vec(i) = first + i - 1
end do
end function seq_unit_stride

pure function cbind_vec_vec(x,y) result(xy)
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: xy(:,:)
integer :: n
n = size(x,1)
if (size(y) /= n) error stop "mismatched sizes in cbind"
xy = reshape([x, y], [n, 2])
end function cbind_vec_vec

pure function cbind_mat_vec(x,y) result(xy)
real(kind=dp), intent(in) :: x(:,:), y(:)
real(kind=dp), allocatable :: xy(:,:)
integer :: n1, n2
n1 = size(x,1)
if (size(y) /= n1) error stop "mismatched sizes in cbind"
n2 = size(x,2)
allocate (xy(n1, n2+1))
xy(:,:n2)  = x
xy(:,n2+1) = y 
end function cbind_mat_vec

pure function cbind_mat_mat(x,y) result(xy)
real(kind=dp), intent(in) :: x(:,:), y(:,:)
real(kind=dp), allocatable :: xy(:,:)
integer :: n1, n2
n1 = size(x,1)
if (size(y,1) /= n1) error stop "mismatched sizes in cbind"
n2 = size(x,2)
allocate (xy(n1, n2+size(y,2)))
xy(:,:n2)  = x
xy(:,n2+1:) = y 
end function cbind_mat_mat

end module util_mod
