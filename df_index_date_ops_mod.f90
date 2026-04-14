module df_index_date_ops_mod
use date_mod, only: date, operator(>=), operator(<=), operator(==), &
   operator(<), operator(>)
use util_mod, only: default
implicit none
private
public :: findloc_index, argsort_index, union_index, intersect_index, &
   is_sorted_index_array, is_unique_index_array, bsearch_exact_index, &
   bsearch_ffill_index, bsearch_bfill_index
contains

pure integer function findloc_index(a, x) result(pos)
! return first location of x in a, or 0 if not found
type(date), intent(in) :: a(:)
type(date), intent(in) :: x
integer :: i
pos = 0
do i = 1, size(a)
   if (a(i) == x) then
      pos = i
      return
   end if
end do
end function findloc_index

subroutine argsort_index(a, perm, ascending)
! return permutation perm such that a(perm) is sorted
 type(date), intent(in) :: a(:)
 integer, allocatable, intent(out) :: perm(:)
 logical, intent(in), optional :: ascending
 logical :: asc
 integer :: n, width, i, left, mid, right, p, q, k
 integer, allocatable :: tmp(:)
 asc = default(.true., ascending)
 n = size(a)
 allocate(perm(n), tmp(n))
 perm = [(i, i=1,n)]
 width = 1
 do while (width < n)
    i = 1
    do while (i <= n)
       left = i
       mid = min(i + width - 1, n)
       right = min(i + 2*width - 1, n)
       p = left
       q = mid + 1
       k = left
       do while (p <= mid .and. q <= right)
          if (asc) then
             if (a(perm(p)) <= a(perm(q))) then
                tmp(k) = perm(p)
                p = p + 1
             else
                tmp(k) = perm(q)
                q = q + 1
             end if
          else
             if (a(perm(p)) >= a(perm(q))) then
                tmp(k) = perm(p)
                p = p + 1
             else
                tmp(k) = perm(q)
                q = q + 1
             end if
          end if
          k = k + 1
       end do
       do while (p <= mid)
          tmp(k) = perm(p)
          p = p + 1
          k = k + 1
       end do
       do while (q <= right)
          tmp(k) = perm(q)
          q = q + 1
          k = k + 1
       end do
       perm(left:right) = tmp(left:right)
       i = i + 2*width
    end do
    width = 2*width
 end do
 deallocate(tmp)
end subroutine argsort_index

pure logical function is_sorted_index_array(a, ascending) result(is_sorted)
! return true if a is sorted
 type(date), intent(in) :: a(:)
 logical, intent(in), optional :: ascending
 logical :: asc
 integer :: i
 asc = default(.true., ascending)
 is_sorted = .true.
 if (size(a) <= 1) return
 if (asc) then
    do i = 2, size(a)
       if (a(i) < a(i-1)) then
          is_sorted = .false.
          return
       end if
    end do
 else
    do i = 2, size(a)
       if (a(i) > a(i-1)) then
          is_sorted = .false.
          return
       end if
    end do
 end if
end function is_sorted_index_array

logical function is_unique_index_array(a) result(is_unique)
! return true if a has no duplicates
 type(date), intent(in) :: a(:)
 integer :: i, n
 integer, allocatable :: perm(:)
 type(date), allocatable :: tmp(:)
 n = size(a)
 is_unique = .true.
 if (n <= 1) return
 allocate(tmp(n))
 tmp = a
 call argsort_index(tmp, perm, ascending=.true.)
 tmp = tmp(perm)
 do i = 2, n
    if (tmp(i) == tmp(i-1)) then
       is_unique = .false.
       exit
    end if
 end do
 deallocate(tmp, perm)
end function is_unique_index_array

function union_index(a, b) result(c)
! return union of a and b preserving first appearance order
 type(date), intent(in) :: a(:), b(:)
 type(date), allocatable :: c(:)
 type(date), allocatable :: tmp(:)
 integer :: n, i
 allocate(tmp(size(a) + size(b)))
 n = 0
 do i = 1, size(a)
    n = n + 1
    tmp(n) = a(i)
 end do
 do i = 1, size(b)
    if (.not. any(tmp(1:n) == b(i))) then
       n = n + 1
       tmp(n) = b(i)
    end if
 end do
 allocate(c(n))
 c = tmp(1:n)
end function union_index

function intersect_index(a, b) result(c)
! return intersection of a and b preserving order from a
 type(date), intent(in) :: a(:), b(:)
 type(date), allocatable :: c(:)
 type(date), allocatable :: tmp(:)
 integer :: n, i
 allocate(tmp(size(a)))
 n = 0
 do i = 1, size(a)
    if (any(b == a(i))) then
      n = n + 1
      tmp(n) = a(i)
    end if
 end do
 allocate(c(n))
 c = tmp(1:n)
end function intersect_index

pure integer function bsearch_exact_index(a, x) result(pos)
! return exact match position in sorted ascending array
 type(date), intent(in) :: a(:)
 type(date), intent(in) :: x
 integer :: lo, hi, mid
 pos = 0
 lo = 1
 hi = size(a)
 do while (lo <= hi)
    mid = (lo + hi)/2
    if (a(mid) == x) then
       pos = mid
       return
    else if (a(mid) < x) then
       lo = mid + 1
    else
       hi = mid - 1
    end if
 end do
end function bsearch_exact_index

pure integer function bsearch_ffill_index(a, x) result(pos)
! return rightmost a(pos) <= x in sorted ascending array
 type(date), intent(in) :: a(:)
 type(date), intent(in) :: x
 integer :: lo, hi, mid
 pos = 0
 lo = 1
 hi = size(a)
 do while (lo <= hi)
    mid = (lo + hi)/2
    if (a(mid) <= x) then
       pos = mid
       lo = mid + 1
    else
       hi = mid - 1
    end if
 end do
end function bsearch_ffill_index

pure integer function bsearch_bfill_index(a, x) result(pos)
! return leftmost a(pos) >= x in sorted ascending array
 type(date), intent(in) :: a(:)
 type(date), intent(in) :: x
 integer :: lo, hi, mid
 pos = 0
 lo = 1
 hi = size(a)
 do while (lo <= hi)
    mid = (lo + hi)/2
    if (a(mid) >= x) then
       pos = mid
       hi = mid - 1
    else
       lo = mid + 1
    end if
 end do
end function bsearch_bfill_index

end module df_index_date_ops_mod
