MODULE inv_cx2
!..use and access
use inv_dp
use inv_core
use inv_mg2
use inv_cxx
implicit none

private
!..procedures
public :: cx_b2, cx_f2, cx_f02

!..data
integer, parameter, public :: &
  cx_nb2(-1:ubound(mg2_nb,dim=1))=(/0,mg2_nb(0:)/)

!..procedures
CONTAINS

SUBROUTINE cx_b2(nki, ik, pc, r, w)
integer, intent (in) :: nki(0:), ik(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=Dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
integer :: i0, i1
real (kind=Dp) :: t0, r0(0:mg2_nk-1,0:mg2_nk-1), &
  y0(0:mg2_nk-1,0:mg2_nk-1), &
  w0(0:cx_dim(mg2_nb,pc%dg)-1), w1(0:cx_dim(mg2_nb,pc%dg)-1)

if (any(size(nki).lt.ik)) then
    stop 'cx_b2: bad nki, ik'
else if (size(ik).ne.mg2_nkk) then
    stop 'cx_b2: bad dimension ik'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
    stop 'cx_b2: bad dimension nki, r'
else if (size(w).ne.cx_dim(mg2_nb,pc%dg)) then
    stop 'cx_b2: bad dimension w'
endif

w1 = 0
if (0.le.pc%dg) then
    do i1 = sum(nki(0:ik(0)-1))+1, sum(nki(0:ik(0)))-1
        do i0 = sum(nki(0:ik(0)-1)), i1-1
            r0 = r((/i0,i1/),(/i0,i1/))
            t0 = cx_cut(pc,r0)
            if (dabs(t0).gt.1d-12) then
                call cx_var(pc, r0, y0)
                call mg2_base(pc%dg, y0, w0)
                w1 = w1 + w0 * t0
            endif
        enddo
    enddo
endif
w = w1

return
END SUBROUTINE cx_b2

FUNCTION cx_f2(nki, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:), cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b2 (nki, (/0/), pc, r, w)
f = dot_product(cf,w)

return
END FUNCTION cx_f2

FUNCTION cx_f02(nki, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:), cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b2 (nki, (/1/), pc, r, w)
f = dot_product(cf,w)

return
END FUNCTION cx_f02

END MODULE inv_cx2
