module kind_mod
    use iso_fortran_env, only: int64
    implicit none
    public :: dp, long_int
    integer, parameter :: dp = kind(1.0d0), long_int=int64
end module kind_mod
