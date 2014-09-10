MODULE inv_mg21
!..use and access
use inv_core
use inv_mgx
implicit none

private
public :: mg21_prims, mg21_prib, mg21_secs, mg21_base
!..data
integer, parameter, private :: nkk=2, nk=3, nr=nk*(nk-1)/2, &
  nkj(0:nkk-1)=(/2,1/)
integer, parameter, public :: &
  mg21_id=5, &
  mg21_nkk=nkk, mg21_nk=nk, mg21_nr=nr, mg21_ngrp=2, &
  mg21_nkj(0:nkk-1)=nkj, &
  mg21_dnpr(0:19) = (/ 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
  mg21_npr(0:19) = (/ 0, 2, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /), &
  mg21_dnpb(0:19) = (/ 1, 2, 4, 6, 9, 12, 16, 20, 25, 30, &
    36, 42, 49, 56, 64, 72, 81, 90, 100, 110 /), &
  mg21_npb(0:19) = (/ 1, 3, 7, 13, 22, 34, 50, 70, 95, 125, &
    161, 203, 252, 308, 372, 444, 525, 615, 715, 825 /), &
  mg21_dnsc(0:19) = (/ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
  mg21_nsc(0:19) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), &
  mg21_dnb(0:19) = (/ 1, 2, 4, 6, 9, 12, 16, 20, 25, 30, &
    36, 42, 49, 56, 64, 72, 81, 90, 100, 110 /), &
  mg21_nb(0:19) = (/ 1, 3, 7, 13, 22, 34, 50, 70, 95, 125, &
    161, 203, 252, 308, 372, 444, 525, 615, 715, 825 /)

!..procedures
CONTAINS

SUBROUTINE mg21_prims(r, u)
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------
integer, parameter :: m=2, m2=m*(m-1)/2
real (kind=dp) :: x(0:nr-1), u0(0:m2-1), u1(0:m-1)

if (size(r,1).ne.nk.or.size(r,2).ne.nk.or.size(u).ne.nr) then
    stop 'mg21_prims: bad dimensions'
endif

!d(i,j)-> x(k)
call mgx_mk1d(nkj, r, x)
call mgx_mk1d(nkj, r, x)
call cg2_prims(x(0:m2-1), u0)
call sym2_prims(x(m2:m2+m-1), u1)
u = (/ u0(0), u1(0), u1(1) /)

return
END SUBROUTINE mg21_prims

SUBROUTINE mg21_prib(mxd, u, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: u(0:)
real (kind=dp), intent (out) :: w(0:mg21_npb(mxd)-1)
!-----------------------------------------------------------------------

if (size(u).ne.nr) then
    stop 'mg21_prib: bad size u'
!! else if (size(w).ne.mg21_npb(mxd)) then
!!  stop 'mg21_prib: bad size w'
endif

! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 2 1 0 0 0 0 0 0 0 &
!   0 0 0 0 0 0 0 0 0 0
! prib,  dnpb(0:*): 1 2 4 6 9 12 16 20 25 30 &
!   36 42 49 56 64 72 81 90 100 110
! constant term
w(0) = 1
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
! terms of degree 10
if (10.le.mxd) then
 w(125:154) = u(0)*w(95:124)
 w(155:159) = u(1)*w(120:124)
 w(160) = u(2)*w(94)
endif
! terms of degree 11
if (11.le.mxd) then
 w(161:196) = u(0)*w(125:160)
 w(197:202) = u(1)*w(155:160)
endif
! terms of degree 12
if (12.le.mxd) then
 w(203:244) = u(0)*w(161:202)
 w(245:250) = u(1)*w(197:202)
 w(251) = u(2)*w(160)
endif
! terms of degree 13
if (13.le.mxd) then
 w(252:300) = u(0)*w(203:251)
 w(301:307) = u(1)*w(245:251)
endif
! terms of degree 14
if (14.le.mxd) then
 w(308:363) = u(0)*w(252:307)
 w(364:370) = u(1)*w(301:307)
 w(371) = u(2)*w(251)
endif
! terms of degree 15
if (15.le.mxd) then
 w(372:435) = u(0)*w(308:371)
 w(436:443) = u(1)*w(364:371)
endif
! terms of degree 16
if (16.le.mxd) then
 w(444:515) = u(0)*w(372:443)
 w(516:523) = u(1)*w(436:443)
 w(524) = u(2)*w(371)
endif
! terms of degree 17
if (17.le.mxd) then
 w(525:605) = u(0)*w(444:524)
 w(606:614) = u(1)*w(516:524)
endif
! terms of degree 18
if (18.le.mxd) then
 w(615:704) = u(0)*w(525:614)
 w(705:713) = u(1)*w(606:614)
 w(714) = u(2)*w(524)
endif
! terms of degree 19
if (19.le.mxd) then
 w(715:814) = u(0)*w(615:714)
 w(815:824) = u(1)*w(705:714)
endif

END SUBROUTINE mg21_prib

SUBROUTINE mg21_secs(mxd, r, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: v(0:)
!-----------------------------------------------------------------------

if (size(r).ne.nk*nk.or.size(v).ne.mg21_nsc(mxd)) then
    stop 'mg21_secs: bad dimensions'
endif

v(0) = 1
! There are no further secondaries

return
END SUBROUTINE mg21_secs

SUBROUTINE mg21_base(mxd, r, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:mg21_nb(mxd)-1)
!-----------------------------------------------------------------------
integer :: k, l0, m0, l, n, ind
real (kind=dp) :: u(0:nr-1), v(0:mg21_nsc(mxd)-1), &
  w0(0:mg21_npb(mxd)-1)

if (size(r,1).ne.nk.or.size(r,2).ne.nk) then
    stop 'mg21_base: bad size r'
!! else if (size(w).ne.mg21_nb(mxd)) then
!!  stop 'mg21_base: bad size w'
endif

call mg21_prims(r, u)
call mg21_prib(mxd, u, w0)
call mg21_secs(mxd, r, v)
! trivial code for this special case
if (size(v).eq.1) then
    w = v(0)*w0
else
    stop 'mg21_base: bad element count'
endif

return
END SUBROUTINE mg21_base

END MODULE inv_mg21
