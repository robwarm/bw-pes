MODULE inv_cx1
!..use and access
use inv_core
use inv_mg1
use inv_cxx
implicit none

private
!..procedures
public :: cx_b1, cx_f1, cx_f01

!..data
integer, parameter, public :: &
  cx_nb1=1

!..procedures
CONTAINS

SUBROUTINE cx_b1(nki, ik, w)
integer, intent (in) :: nki(0:), ik(0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------

if (any(size(nki).lt.ik)) then
    stop 'cx_b1: bad nki, ik'
else if (size(ik).ne.mg1_nkk) then
    stop 'cx_b1: bad dimension ik'
else if (size(w).ne.1) then
    stop 'cx_b1: bad dimension w'
endif

w = (/ nki(ik(0)) /)

END SUBROUTINE cx_b1

FUNCTION cx_f1(nki, cf) RESULT (f)
integer, intent (in) :: nki(0:)
real (kind=dp), intent (in) :: cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b1(nki, (/0/), w)
f = dot_product(cf,w)

return
END FUNCTION cx_f1

FUNCTION cx_f01(nki, cf) RESULT (f)
integer, intent (in) :: nki(0:)
real (kind=dp), intent (in) :: cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b1(nki,(/1/), w)
f = dot_product(cf,w)

return
END FUNCTION cx_f01

END MODULE inv_cx1
