MODULE inv_mg11
!..use and access
use inv_dp
use inv_core
use inv_mgx
implicit none

private
public :: mg11_prims, mg11_prib, mg11_secs, mg11_base

!..data
! Note: all rather trivial.  We keep it this way to be consistent
! with other routines in the inv package.
integer, private :: i_local
integer, parameter, private :: nkk=2, nk=2, nr=nk*(nk-1)/2, &
  nkj(0:nkk-1)=(/1,1/)
integer, parameter, public :: &
  mg11_id=3, &
  mg11_nkk=nkk, mg11_nk=nk, mg11_nr=nr, mg11_ngrp=1, &
  mg11_nkj(0:nkk-1)=nkj, &
  mg11_dnpr(0:19) = (/ 0, 1, (0,i_local=2,19) /), &
  mg11_npr(0:19) = (/ 0, (1,i_local=1,19) /), &
  mg11_dnpb(0:19) = (/ (1,i_local=0,19) /), &
  mg11_npb(0:19) = (/ (i_local+1,i_local=0,19) /), &
  mg11_dnsc(0:19) = (/ 1, (0,i_local=1,19) /), &
  mg11_nsc(0:19) = (/ (1,i_local=0,19) /), &
  mg11_dnb(0:19) = (/ (1,i_local=0,19) /), &
  mg11_nb(0:19) = (/ (i_local+1,i_local=0,19) /)

!..procedures
CONTAINS

SUBROUTINE mg11_prims(r, u)
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------

if (size(r,1).ne.nk.or.size(r,2).ne.nk.or.size(u).ne.nr) then
    stop 'mg11_prims: bad dimensions'
endif

! There is just one variable
u(0) = (r(0,1)+r(1,0))/2

return
END SUBROUTINE mg11_prims

SUBROUTINE mg11_prib(mxd, u, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: u(0:)
real (kind=dp), intent (out) :: w(0:mg11_npb(mxd)-1)
!-----------------------------------------------------------------------
integer :: i

if (size(u).ne.nr) then
    stop 'mg11_prib: bad size u'
!! else if (size(w).ne.mg11_npb(mxd)) then
!!  stop 'mg11_prib: bad size w'
endif

! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 1 0 0 0 0 0 0 0 0
! prib,  dnpb(0:*): 1 1 1 1 1 1 1 1 1 1
! constant term
w(0) = 1
do i=1, mxd
    w(i) = u(0) * w(i-1)
enddo

END SUBROUTINE mg11_prib

SUBROUTINE mg11_secs(mxd, r, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
!! real (kind=dp), intent (out) :: v(0:)
real (kind=dp), intent (out) :: v(0:0)
!-----------------------------------------------------------------------

if (size(r).ne.nk*nk.or.size(v).ne.mg11_nsc(mxd)) then
    stop 'mg11_secs: bad dimensions'
endif

! There is only the trivial secondary
v(0) = 1
return

END SUBROUTINE mg11_secs

SUBROUTINE mg11_base(mxd, r, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:mg11_nb(mxd)-1)
!-----------------------------------------------------------------------
integer :: k, l0, m0, l, n, ind
real (kind=dp) :: u(0:nr-1), v(0:mg11_nsc(mxd)-1), &
  w0(0:mg11_npb(mxd)-1)

if (size(r,1).ne.nk.or.size(r,2).ne.nk) then
    stop 'mg11_base: bad size r'
!! else if (size(w).ne.mg11_nb(mxd)) then
!!  stop 'mg11_base: bad size w'
endif

call mg11_prims(r, u)
call mg11_prib(mxd, u, w0)
call mg11_secs(mxd, r, v)
! trivial code for this special case
if (size(v).eq.1) then
    w = v(0)*w0
else
    stop 'mg11_base: bad element count'
endif

return
END SUBROUTINE mg11_base

END MODULE inv_mg11
