MODULE inv_cx21
!..use and access
use inv_core
use inv_mg21
use inv_cxx
implicit none

private
!..procedures
public :: cx_b21, cx_f21, cx_f12, cx_b21_consec, cx_f21_consec

!..data
integer, parameter, public :: &
  cx_nb21(-1:ubound(mg21_nb,dim=1))=(/0,mg21_nb(0:)/)

!..procedures
CONTAINS

SUBROUTINE cx_b21(nki, ik, pc, r, w)
integer, intent (in) :: nki(0:), ik(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
integer :: i0, i1, j0
real (kind=dp) :: t0, r0(0:mg21_nk-1,0:mg21_nk-1), &
  y0(0:mg21_nk-1,0:mg21_nk-1), &
  w0(0:cx_dim(mg21_nb,pc%dg)-1), w1(0:cx_dim(mg21_nb,pc%dg)-1)

if (any(size(nki).lt.ik)) then
    stop 'cx_b21: bad nki, ik'
else if (size(ik).ne.mg21_nkk) then
    stop 'cx_b21: bad dimension ik'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
    stop 'cx_b21: bad dimension nki, r'
else if (size(w).ne.cx_dim(mg21_nb,pc%dg)) then
    stop 'cx_b21: bad dimension w'
endif

w1 = 0
if (0.le.pc%dg) then
    do j0 = sum(nki(0:ik(1)-1)), sum(nki(0:ik(1)))-1
        do i1 = sum(nki(0:ik(0)-1))+1, sum(nki(0:ik(0)))-1
            do i0 = sum(nki(0:ik(0)-1)), i1-1
                r0 = r((/i0,i1,j0/),(/i0,i1,j0/))
                t0 = cx_cut(pc,r0)
                if (t0.ne.0) then
                    call cx_var(pc, r0, y0)
                    call mg21_base(pc%dg, y0, w0)
                    w1 = w1 + w0 * t0
                endif
            enddo
        enddo
    enddo
endif
w = w1

return
END SUBROUTINE cx_b21

SUBROUTINE cx_b21_consec(nki, ik, pc, cp, r, w, whichlin)
integer, intent (in) :: nki(0:), ik(0:)
type (cx_t_consec), intent (in) :: pc
type (cp_t_consec), dimension(2), intent(in) :: cp
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:)
integer, intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: i0, i1, j0
real (kind=dp) :: t0, r0(0:mg21_nk-1,0:mg21_nk-1), &
  y0(0:mg21_nk-1,0:mg21_nk-1), &
  w0(0:cx_dim(mg21_nb,pc%dg)-1), w1(0:cx_dim(mg21_nb,pc%dg)-1)

if (any(size(nki).lt.ik)) then
    stop 'cx_b21_consec: bad nki, ik'
else if (size(ik).ne.mg21_nkk) then
    stop 'cx_b21_consec: bad dimension ik'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
    stop 'cx_b21_consec: bad dimension nki, r'
else if (size(w).ne.cx_dim(mg21_nb,pc%dg)) then
    stop 'cx_b21_consec: bad dimension w'
endif

w1 = 0
if (0.le.pc%dg) then
    do j0 = sum(nki(0:ik(1)-1)), sum(nki(0:ik(1)))-1
        do i1 = sum(nki(0:ik(0)-1))+1, sum(nki(0:ik(0)))-1
            do i0 = sum(nki(0:ik(0)-1)), i1-1
                r0 = r((/i0,i1,j0/),(/i0,i1,j0/))
                t0 = cx_cut_consec(pc,r0,whichlin)
                if (t0.ne.0) then
                    call cx_var_consec(pc, cp, r0, y0, whichlin)
                    call mg21_base(pc%dg, y0, w0)
                    w1 = w1 + w0 * t0
                endif
            enddo
        enddo
    enddo
endif
w = w1

return
END SUBROUTINE cx_b21_consec

FUNCTION cx_f21(nki, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:), cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b21(nki, (/0,1/), pc, r, w)
f = dot_product(cf,w)

return
END FUNCTION cx_f21

FUNCTION cx_f21_consec(nki, r, pc, cf, cp, whichlin) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t_consec), intent (in) :: pc
type (cp_t_consec), dimension(2), intent(in) :: cp
real (kind=dp), intent (in) :: r(0:,0:), cf(0:)
real (kind=dp) :: f
integer, intent(in) :: whichlin
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b21_consec (nki, (/0,1/), pc, cp, r, w,whichlin)
f = dot_product(cf,w)

return
END FUNCTION cx_f21_consec

FUNCTION cx_f12(nki, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:), cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_b21 (nki, (/1,0/), pc, r, w)
f = dot_product(cf,w)

return
END FUNCTION cx_f12

END MODULE inv_cx21
