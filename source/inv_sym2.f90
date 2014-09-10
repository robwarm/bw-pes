MODULE inv_sym2
!..use and access
use inv_dp
implicit none

private
public :: sym2_gens, sym2_prims, sym2_secs, sym2_base

!..data
integer, parameter, public :: &
  sym2_nr=2, sym2_ngrp=2, sym2_ngen=1, &
  sym2_dnpr(0:9) = (/ 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 /), &
  sym2_npr(0:9) = (/ 0, 1, 2, 2, 2, 2, 2, 2, 2, 2 /), &
  sym2_dnsc(0:9) = (/ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
  sym2_nsc(0:9) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), &
  sym2_dnb(0:9) = (/ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 /), &
  sym2_nb(0:9) = (/ 1, 2, 4, 6, 9, 12, 16, 20, 25, 30 /)

!..procedures
CONTAINS

SUBROUTINE sym2_gens(ind, iord)
integer, intent (in) :: ind
integer, intent (out) :: iord(0:)
!-----------------------------------------------------------------------

if (size(iord).ne.sym2_nr) then
    stop 'sym2_gens: bad size iord'
endif

select case (ind)
    case (0)
        ! permutation (0,1)
        iord = (/ 1, 0 /)
    case default
        stop 'sym2_gens: invalid index'
end select

return
END SUBROUTINE sym2_gens

SUBROUTINE sym2_prims(x, u)
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------

if (size(x).ne.sym2_nr.or.size(u).ne.sym2_nr) then
    stop 'sym2_prims: bad dimensions'
endif

u(0) = sum(x) / size(x)
u(1) = sum(x**2) / size(x)

return
END SUBROUTINE sym2_prims

SUBROUTINE sym2_secs(mxd, x, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: v(0:)
!-----------------------------------------------------------------------

if (size(x).ne.sym2_nr.or.size(v).ne.sym2_nsc(mxd)) then
    stop 'sym2_secs: bad dimensions'
endif

! There is only the trivial secondary
v(0) = 1

return
END SUBROUTINE sym2_secs

SUBROUTINE sym2_base(mxd, x, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
real (kind=dp) :: u(0:sym2_nr-1), v(0:sym2_nsc(mxd)-1)

if (size(x).ne.sym2_nr.or.size(w).ne.sym2_nb(mxd)) then
    stop 'sym2_base: bad dimensions'
endif

call sym2_prims(x, u)
call sym2_secs(mxd, x, v)
! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 1 1 0 0 0 0 0 0 0
! secs,  dnsc(0:*): 1 0 0 0 0 0 0 0 0 0
! base,  dnb(0:*): 1 1 2 2 3 3 4 4 5 5
! constant term
w(0) = v(0)
! terms of degree 1
if (1.le.mxd) then
    w(1) = u(0)*w(0)
endif
! terms of degree 2
if (2.le.mxd) then
    w(2) = u(0)*w(1)
    w(3) = u(1)*w(0)
endif
! terms of degree 3
if (3.le.mxd) then
    w(4:5) = u(0)*w(2:3)
endif
! terms of degree 4
if (4.le.mxd) then
    w(6:7) = u(0)*w(4:5)
    w(8) = u(1)*w(3)
endif
! terms of degree 5
if (5.le.mxd) then
    w(9:11) = u(0)*w(6:8)
endif
! terms of degree 6
if (6.le.mxd) then
    w(12:14) = u(0)*w(9:11)
    w(15) = u(1)*w(8)
endif
! terms of degree 7
if (7.le.mxd) then
    w(16:19) = u(0)*w(12:15)
endif
! terms of degree 8
if (8.le.mxd) then
    w(20:23) = u(0)*w(16:19)
    w(24) = u(1)*w(15)
endif
! terms of degree 9
if (9.le.mxd) then
    w(25:29) = u(0)*w(20:24)
endif

END SUBROUTINE sym2_base

END MODULE inv_sym2
