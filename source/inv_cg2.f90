MODULE inv_cg2
!..use and access
use inv_dp
implicit none

private
public :: cg2_prims, cg2_prib, cg2_secs, cg2_base

!..data
! Note: all rather trivial.  We keep it this way to be consistent
! with other routines in the inv package.
integer, private :: i_local
integer, parameter, private :: nk=2, nr=nk*(nk-1)/2
integer, parameter, public :: &
  cg2_nk=nk, cg2_nr=nr, cg2_ngrp=2, cg2_ngen=1, &
  cg2_dnpr(0:19) = (/ 0, 1, (0,i_local=2,19) /), &
  cg2_npr(0:19) = (/ 0, (1,i_local=1,19) /), &
  cg2_dnpb(0:19) = (/ (1,i_local=0,19) /), &
  cg2_npb(0:19) = (/ (i_local+1,i_local=0,19) /), &
  cg2_dnsc(0:19) = (/ 1, (0,i_local=1,19) /), &
  cg2_nsc(0:19) = (/ (1,i_local=0,19) /), &
  cg2_dnb(0:19) = (/ (1,i_local=0,19) /), &
  cg2_nb(0:19) = (/ (i_local+1,i_local=0,19) /)

!..procedures
interface cg2_base
    module procedure cg2_base_vec, cg2_base_scal
end interface

CONTAINS

SUBROUTINE cg2_prims(x, u)
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: u(0:)
!-----------------------------------------------------------------------

if (size(x).ne.nr.or.size(u).ne.nr) then
    stop 'cg2_prims: bad dimensions'
endif

! There is just one variable
u = x

return
END SUBROUTINE cg2_prims

SUBROUTINE cg2_prib(mxd, x, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: w(0:cg2_npb(mxd)-1)
!-----------------------------------------------------------------------
real (kind=dp) :: u(0:nr-1)
integer :: i
if (size(x).ne.nr) then
    stop 'cg2_prib: bad size x'
!! else if (size(w).ne.cg2_npb(mxd)) then
!!  stop 'cg2_prib: bad size w'
endif

call cg2_prims(x, u)
! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 1 0 0 0 0 0 0 0 0
! prib,  dnpb(0:*): 1 1 1 1 1 1 1 1 1 1
! empty product
w(0) = 1
do i=1, mxd
    w(i) = u(0)*w(i-1)
enddo

END SUBROUTINE cg2_prib

SUBROUTINE cg2_secs(mxd, x, v)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:)
!! real (kind=dp), intent (out) :: v(0:)
real (kind=dp), intent (out) :: v(0:0)
!-----------------------------------------------------------------------

if (size(x).ne.nr.or.size(v).ne.cg2_nsc(mxd)) then
    stop 'cg2_secs: bad dimensions'
endif

! There is only the trivial secondary
v(0) = 1

return
END SUBROUTINE cg2_secs

SUBROUTINE cg2_base_vec(mxd, x, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x(0:nr-1)
real (kind=dp), intent (out) :: w(0:1+mxd-1)
!-----------------------------------------------------------------------
real (kind=dp) :: u(0:nr-1), v(0:0)
integer :: i

call cg2_prims(x, u)
call cg2_secs(mxd, x, v)
! Trivial code.  Different from other routines in the inv package,
! this one does not have an upper restriction on mxd.  There is only
! one primary and only one secondary.
! constant term
w(0) = v(0)
! higher degrees
do i = 1, mxd
    w(i) = u(0)*w(i-1)
enddo

END SUBROUTINE cg2_base_vec

SUBROUTINE cg2_base_scal(mxd, x0, w)
integer, intent (in) :: mxd
real (kind=dp), intent (in) :: x0
real (kind=dp), intent (out) :: w(0:1+mxd-1)
!-----------------------------------------------------------------------
integer :: i
! Trivial code.  Different from other routines in the inv package,
! this one does not have an upper restriction on mxd.
! constant term
w(0) = 1
! higher degrees
do i = 1, mxd
    w(i) = x0*w(i-1)
enddo

END SUBROUTINE cg2_base_scal

END MODULE inv_cg2
