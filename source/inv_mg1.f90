MODULE inv_mg1
!..use and access
use inv_dp
use inv_core
use inv_mgx
implicit none

private
public :: mg1_prims, mg1_prib, mg1_secs, mg1_base

!..data
! Note: all totally trivial.  We keep it this way to be consistent
! with other routines in the inv package.
integer, private :: i_local
integer, parameter, private :: nkk=1, nk=1, nr=nk*(nk-1)/2, &
  nkj(0:nkk-1)=(/1/)
integer, parameter, public :: &
  mg1_id=1, &
  mg1_nkk=nkk, mg1_nk=nk, mg1_nr=nr, mg1_ngrp=1, &
  mg1_nkj(0:nkk-1)=nkj, &
  mg1_dnpr(0:19) = 0, &
  mg1_npr(0:19) = 0, &
  mg1_dnpb(0:19) = (/ 1, (0,i_local=1,19) /), &
  mg1_npb(0:19) = 1, &
  mg1_dnsc(0:19) = (/ 1, (0,i_local=1,19) /), &
  mg1_nsc(0:19) = 1, &
  mg1_dnb(0:19) = (/ 1, (0,i_local=1,19) /), &
  mg1_nb(0:19) = 1

  !..procedures
CONTAINS

SUBROUTINE mg1_prims(r, u)
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------
integer :: i

if (size(r,1).ne.nk.or.size(r,2).ne.nk.or.size(u).ne.nr) then
    stop 'mg1_prims: bad dimensions'
endif

! There are no primaries
u = (/ (0.0_dp,i=1,0) /) ! an empty array of well-defined type

return
END SUBROUTINE mg1_prims

SUBROUTINE mg1_prib(mxd, u, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: u(0:)
real (kind=dp), intent (out) :: w(0:mg1_npb(mxd)-1)
!-----------------------------------------------------------------------

if (size(u).ne.nr) then
    stop 'mg1_prib: bad size u'
!! else if (size(w).ne.mg1_npb(mxd)) then
!!  stop 'mg1_prib: bad size w'
endif

! constant term
w(0) = 1

END SUBROUTINE mg1_prib

SUBROUTINE mg1_secs(mxd, r, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
!! real (kind=dp), intent (out) :: v(0:)
real (kind=dp), intent (out) :: v(0:0)
!-----------------------------------------------------------------------

if (size(r).ne.nk*nk.or.size(v).ne.mg1_nsc(mxd)) then
    stop 'mg1_secs: bad dimensions'
endif

! There is only the trivial secondary
v(0) = 1

return
END SUBROUTINE mg1_secs

SUBROUTINE mg1_base(mxd, r, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: r(0:,0:)
!! real (kind=wp), intent (out) :: w(0:mg1_nb(mxd)-1)
real (kind=dp), intent (out) :: w(0:0)
!-----------------------------------------------------------------------

if (size(r,1).ne.nk.or.size(r,2).ne.nk) then
    stop 'mg1_base: bad size r'
!! else if (size(w).ne.mg1_nb(mxd)) then
!!  stop 'mg1_base: bad size w'
endif

! trivial code for this special case
call mg1_secs(mxd, r, w)

return
END SUBROUTINE mg1_base

END MODULE inv_mg1
