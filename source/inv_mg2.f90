MODULE inv_mg2
!..use and access
use inv_core
use inv_mgx
implicit none

private
public :: mg2_prims, mg2_prib, mg2_secs, mg2_base

!..data
! Note: all rather trivial.  We keep it this way to be consistent
! with other routines in the inv package.
integer, private :: i_local
integer, parameter, private :: nkk=1, nk=2, nr=nk*(nk-1)/2, &
  nkj(0:nkk-1)=(/2/)
integer, parameter, public :: &
  mg2_id=2, &
  mg2_nkk=nkk, mg2_nk=nk, mg2_nr=nr, mg2_ngrp=2, &
  mg2_nkj(0:nkk-1)=nkj, &
  mg2_dnpr(0:19) = (/ 0, 1, (0,i_local=2,19) /), &
  mg2_npr(0:19) = (/ 0, (1,i_local=1,19) /), &
  mg2_dnpb(0:19) = (/ (1,i_local=0,19) /), &
  mg2_npb(0:19) = (/ (i_local+1,i_local=0,19) /), &
  mg2_dnsc(0:19) = (/ 1, (0,i_local=1,19) /), &
  mg2_nsc(0:19) = (/ (1,i_local=0,19) /), &
  mg2_dnb(0:19) = (/ (1,i_local=0,19) /), &
  mg2_nb(0:19) = (/ (i_local+1,i_local=0,19) /)

!..procedures
CONTAINS

SUBROUTINE mg2_prims(r, u)
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------

if (size(r,1).ne.nk.or.size(r,2).ne.nk.or.size(u).ne.nr) then
    stop 'mg2_prims: bad dimensions'
endif

! There is just one variable
u(0) = (r(0,1)+r(1,0))/2

return
END SUBROUTINE mg2_prims

SUBROUTINE mg2_prib(mxd, u, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: u(0:)
real (kind=dp), intent (out) :: w(0:mg2_npb(mxd)-1)
!-----------------------------------------------------------------------

if (size(u).ne.nr) then
    stop 'mg2_prib: bad size u'
!! else if (size(w).ne.mg2_npb(mxd)) then
!!  stop 'mg2_prib: bad size w'
endif

! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 1 0 0 0 0 0 0 0 0
! prib,  dnpb(0:*): 1 1 1 1 1 1 1 1 1 1
! constant term
w(0) = 1
do i=1, mxd
    w(i) = u(0) * w(i-1)
enddo

END SUBROUTINE mg2_prib

SUBROUTINE mg2_secs(mxd, r, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
!! real (kind=wp), intent (out) :: v(0:)
real (kind=dp), intent (out) :: v(0:0)
!-----------------------------------------------------------------------

if (size(r).ne.nk*nk.or.size(v).ne.mg2_nsc(mxd)) then
    stop 'mg2_secs: bad dimensions'
endif

! There is only the trivial secondary
v(0) = 1

return
END SUBROUTINE mg2_secs

SUBROUTINE mg2_base(mxd, r, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:mg2_nb(mxd)-1)
!-----------------------------------------------------------------------
integer :: k, l0, m0, l, n, ind
real (kind=wp) :: u(0:nr-1), v(0:mg2_nsc(mxd)-1), w0(0:mg2_npb(mxd)-1)

if (size(r,1).ne.nk.or.size(r,2).ne.nk) then
    stop 'mg2_base: bad size r'
!! else if (size(w).ne.mg2_nb(mxd)) then
!!  stop 'mg2_base: bad size w'
endif

call mg2_prims(r, u)
call mg2_prib(mxd, u, w0)
call mg2_secs(mxd, r, v)
! trivial code for this special case
if (size(v).eq.1) then
    w = v(0)*w0
else
    stop 'mg2_base: bad element count'
endif

return
END SUBROUTINE mg2_base


END MODULE inv_mg2
