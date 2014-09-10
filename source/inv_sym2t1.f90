MODULE inv_sym2t1
!..use and access
use inv_dp
implicit none

private
public :: sym2t1_gens, sym2t1_prims, sym2t1_secs, sym2t1_base

!..data
integer, parameter, public :: &
  sym2t1_nr=3, sym2t1_ngrp=2, sym2t1_ngen=1, &
  sym2t1_dnpr(0:9) = (/ 0, 2, 1, 0, 0, 0, 0, 0, 0, 0 /), &
  sym2t1_npr(0:9) = (/ 0, 2, 3, 3, 3, 3, 3, 3, 3, 3 /), &
  sym2t1_dnsc(0:9) = (/ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
  sym2t1_nsc(0:9) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), &
  sym2t1_dnb(0:9) = (/ 1, 2, 4, 6, 9, 12, 16, 20, 25, 30 /), &
  sym2t1_nb(0:9) = (/ 1, 3, 7, 13, 22, 34, 50, 70, 95, 125 /)

!..procedures
CONTAINS

SUBROUTINE sym2t1_gens(ind, iord)
integer, intent (in) :: ind
integer, intent (out) :: iord(0:)
!-----------------------------------------------------------------------

if (size(iord).ne.sym2t1_nr) then
    stop 'sym2t1_gens: bad size iord'
endif

select case (ind)
    case (0)
        ! permutation (0,1) on the left of the s2 part
        ! (00, 10) =>
        ! (10, 00)
        iord = (/ 1, 0, 2 /)
    case default
        stop 'sym2t1_gens: invalid index'
end select

return
END SUBROUTINE sym2t1_gens

SUBROUTINE sym2t1_prims(x, u)
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------

if (size(x).ne.sym2t1_nr.or.size(u).ne.sym2t1_nr) then
    stop 'sym2t1_prims: bad dimensions'
endif

u(0) = sum(x(0:1))/2
u(1) = x(2)
u(2) = sum(x(0:1)**2)/2

return
END SUBROUTINE sym2t1_prims

SUBROUTINE sym2t1_secs(mxd, x, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: v(0:)
!-----------------------------------------------------------------------

if (size(x).ne.sym2t1_nr.or.size(v).ne.sym2t1_nsc(mxd)) then
    stop 'sym2t1_secs: bad dimensions'
endif

! There is only the trivial secondary
v(0) = 1

return
END SUBROUTINE sym2t1_secs

SUBROUTINE sym2t1_base(mxd, x, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
real (kind=dp) :: u(0:sym2t1_nr-1), v(0:sym2t1_nsc(mxd)-1)

if (size(x).ne.sym2t1_nr.or.size(w).ne.sym2t1_nb(mxd)) then
    stop 'sym2t1_base: bad dimensions'
endif

call sym2t1_prims(x, u)
call sym2t1_secs(mxd, x, v)
! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 2 1 0 0 0 0 0 0 0
! secs,  dnsc(0:*): 1 0 0 0 0 0 0 0 0 0
! base,  dnb(0:*): 1 2 4 6 9 12 16 20 25 30
! constant term
w(0) = v(0)
! terms of degree 1
if (1.le.mxd) then
 w(1) = u(0)*w(0)
 w(2) = u(1)*w(0)
endif
! terms of degree 2
if (2.le.mxd) then
 w(3:4) = u(0)*w(1:2)
 w(5) = u(1)*w(2)
 w(6) = u(2)*w(0)
endif
! terms of degree 3
if (3.le.mxd) then
 w(7:10) = u(0)*w(3:6)
 w(11:12) = u(1)*w(5:6)
endif
! terms of degree 4
if (4.le.mxd) then
 w(13:18) = u(0)*w(7:12)
 w(19:20) = u(1)*w(11:12)
 w(21) = u(2)*w(6)
endif
! terms of degree 5
if (5.le.mxd) then
 w(22:30) = u(0)*w(13:21)
 w(31:33) = u(1)*w(19:21)
endif
! terms of degree 6
if (6.le.mxd) then
 w(34:45) = u(0)*w(22:33)
 w(46:48) = u(1)*w(31:33)
 w(49) = u(2)*w(21)
endif
! terms of degree 7
if (7.le.mxd) then
 w(50:65) = u(0)*w(34:49)
 w(66:69) = u(1)*w(46:49)
endif
! terms of degree 8
if (8.le.mxd) then
 w(70:89) = u(0)*w(50:69)
 w(90:93) = u(1)*w(66:69)
 w(94) = u(2)*w(49)
endif
! terms of degree 9
if (9.le.mxd) then
 w(95:119) = u(0)*w(70:94)
 w(120:124) = u(1)*w(90:94)
endif

END SUBROUTINE sym2t1_base

END MODULE inv_sym2t1
