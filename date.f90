module date_mod
implicit none
private
public :: date, valid, date_from_iso, date_from_basic, &
   operator(+), operator(-), operator(==), operator(/=), &
   operator(<), operator(<=), operator(>), operator(>=)

type :: date
   integer :: year = 0
   integer :: month = 0
   integer :: day = 0
contains
   procedure :: to_str ! return date as yyyy-mm-dd
end type date

interface operator(+)
   module procedure add_days_right
   module procedure add_days_left
end interface

interface operator(-)
   module procedure subtract_days
   module procedure difference_days
end interface

interface operator(==)
   module procedure eq_date
end interface

interface operator(/=)
   module procedure ne_date
end interface

interface operator(<)
   module procedure lt_date
end interface

interface operator(<=)
   module procedure le_date
end interface

interface operator(>)
   module procedure gt_date
end interface

interface operator(>=)
   module procedure ge_date
end interface

contains

pure elemental function to_str(this) result(s) ! return date as yyyy-mm-dd
class(date), intent(in) :: this
character(len=10) :: s
s = zero_pad_4(this%year) // '-' // zero_pad_2(this%month) // '-' // zero_pad_2(this%day)
end function to_str

pure elemental logical function valid(x) ! return true if the date is valid
type(date), intent(in) :: x
valid = .false.
if (x%month < 1 .or. x%month > 12) return
if (x%day < 1) return
if (x%day > days_in_month(x%year, x%month)) return
valid = .true.
end function valid

pure elemental integer function days_in_month(year, month) ! return number of days in a year-month pair
integer, intent(in) :: year, month
select case (month)
case (1,3,5,7,8,10,12)
   days_in_month = 31
case (4,6,9,11)
   days_in_month = 30
case (2)
   if ((mod(year,4) == 0 .and. mod(year,100) /= 0) .or. mod(year,400) == 0) then
      days_in_month = 29
   else
      days_in_month = 28
   end if
case default
   days_in_month = 0
end select
end function days_in_month

pure function date_from_iso(s) result(x) ! convert yyyy-mm-dd to a date
character(len=*), intent(in) :: s
type(date) :: x
character(len=len(s)) :: t
integer :: y, m, d
logical :: ok1, ok2, ok3

x = date(0,0,0)
t = adjustl(s)
if (len_trim(t) /= 10) return
if (t(5:5) /= '-' .or. t(8:8) /= '-') return
call parse_uint(t(1:4), y, ok1)
call parse_uint(t(6:7), m, ok2)
call parse_uint(t(9:10), d, ok3)
if (.not. (ok1 .and. ok2 .and. ok3)) return
x = date(y,m,d)
end function date_from_iso

pure function date_from_basic(s) result(x) ! convert yyyymmdd to a date
character(len=*), intent(in) :: s
type(date) :: x
character(len=len(s)) :: t
integer :: y, m, d
logical :: ok1, ok2, ok3

x = date(0,0,0)
t = adjustl(s)
if (len_trim(t) /= 8) return
call parse_uint(t(1:4), y, ok1)
call parse_uint(t(5:6), m, ok2)
call parse_uint(t(7:8), d, ok3)
if (.not. (ok1 .and. ok2 .and. ok3)) return
x = date(y,m,d)
end function date_from_basic

pure elemental type(date) function add_days_right(x, n) ! return date plus integer number of days
type(date), intent(in) :: x
integer, intent(in) :: n
if (.not. valid(x)) then
   add_days_right = date(0,0,0)
else
   add_days_right = from_day_number(day_number(x) + n)
end if
end function add_days_right

pure elemental type(date) function add_days_left(n, x) ! return integer number of days plus date
integer, intent(in) :: n
type(date), intent(in) :: x
add_days_left = add_days_right(x, n)
end function add_days_left

pure elemental type(date) function subtract_days(x, n) ! return date minus integer number of days
type(date), intent(in) :: x
integer, intent(in) :: n
subtract_days = add_days_right(x, -n)
end function subtract_days

pure elemental integer function difference_days(x, y) ! return number of days between two dates
type(date), intent(in) :: x, y
if (.not. valid(x) .or. .not. valid(y)) then
   difference_days = 0
else
   difference_days = day_number(x) - day_number(y)
end if
end function difference_days

pure elemental logical function eq_date(x, y) ! return true if two dates are equal
type(date), intent(in) :: x, y
eq_date = x%year == y%year .and. x%month == y%month .and. x%day == y%day
end function eq_date

pure elemental logical function ne_date(x, y) ! return true if two dates are not equal
type(date), intent(in) :: x, y
ne_date = .not. eq_date(x, y)
end function ne_date

pure elemental logical function lt_date(x, y) ! return true if left date is earlier than right date
type(date), intent(in) :: x, y
lt_date = x%year < y%year .or. (x%year == y%year .and. (x%month < y%month .or. (x%month == y%month .and. x%day < y%day)))
end function lt_date

pure elemental logical function le_date(x, y) ! return true if left date is earlier than or equal to right date
type(date), intent(in) :: x, y
le_date = lt_date(x, y) .or. eq_date(x, y)
end function le_date

pure elemental logical function gt_date(x, y) ! return true if left date is later than right date
type(date), intent(in) :: x, y
gt_date = .not. le_date(x, y)
end function gt_date

pure elemental logical function ge_date(x, y) ! return true if left date is later than or equal to right date
type(date), intent(in) :: x, y
ge_date = .not. lt_date(x, y)
end function ge_date

pure elemental integer function day_number(x) ! return day count used for arithmetic
type(date), intent(in) :: x
integer :: y, m, d, era, yoe, doy, doe, mp
y = x%year
m = x%month
d = x%day
if (m <= 2) y = y - 1
era = floor_div(y, 400)
yoe = y - era*400
if (m > 2) then
   mp = m - 3
else
   mp = m + 9
end if
doy = (153*mp + 2)/5 + d - 1
doe = yoe*365 + yoe/4 - yoe/100 + doy
day_number = era*146097 + doe - 719468
end function day_number

pure elemental type(date) function from_day_number(z) ! convert internal day count to a date
integer, intent(in) :: z
integer :: zz, era, doe, yoe, y, doy, mp, m, d
zz = z + 719468
era = floor_div(zz, 146097)
doe = zz - era*146097
yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365
y = yoe + era*400
doy = doe - (365*yoe + yoe/4 - yoe/100)
mp = (5*doy + 2)/153
d = doy - (153*mp + 2)/5 + 1
if (mp < 10) then
   m = mp + 3
else
   m = mp - 9
end if
if (m <= 2) y = y + 1
from_day_number = date(y,m,d)
end function from_day_number

pure elemental integer function floor_div(a, b) ! return floor(a/b) for positive b
integer, intent(in) :: a, b
floor_div = a / b
if (mod(a, b) < 0) floor_div = floor_div - 1
end function floor_div

pure elemental function zero_pad_2(n) result(s) ! return a 2-character zero-padded integer string
integer, intent(in) :: n
character(len=2) :: s
if (n < 0 .or. n > 99) then
   s = '**'
   return
end if
s(1:1) = achar(iachar('0') + n/10)
s(2:2) = achar(iachar('0') + mod(n,10))
end function zero_pad_2

pure elemental function zero_pad_4(n) result(s) ! return a 4-character zero-padded integer string
integer, intent(in) :: n
character(len=4) :: s
integer :: m
if (n < 0 .or. n > 9999) then
   s = '****'
   return
end if
m = n
s(4:4) = achar(iachar('0') + mod(m,10))
m = m / 10
s(3:3) = achar(iachar('0') + mod(m,10))
m = m / 10
s(2:2) = achar(iachar('0') + mod(m,10))
m = m / 10
s(1:1) = achar(iachar('0') + mod(m,10))
end function zero_pad_4

pure subroutine parse_uint(s, n, ok) ! parse a nonnegative integer from a string
character(len=*), intent(in) :: s
integer, intent(out) :: n
logical, intent(out) :: ok
integer :: i, m, digit
character(len=len(s)) :: t
n = 0
ok = .false.
t = adjustl(s)
m = len_trim(t)
if (m <= 0) return
do i = 1, m
   if (t(i:i) < '0' .or. t(i:i) > '9') return
   digit = iachar(t(i:i)) - iachar('0')
   n = 10*n + digit
end do
ok = .true.
end subroutine parse_uint

end module date_mod
