MODULE pes_x2y1
!..use and access
use inv
use pes0
use pes_x1
use pes_x2
use pes_y1
use pes_x1y1
implicit none

private
public :: pes_x2y1_read, pes_x2y1_pot, pes_x2y1_pot_r, pes_x2y1_add, pes_x2y1_getcf

!..data
save
integer, parameter, public :: &
  pes_x2y1_nki(0:1)=(/2,1/), pes_x2y1_nk=3, &
  pes_x2y1_nb(-1:size(cx_nb21(-1:))-2)=cx_nb21(-1:)
character (len=*), parameter, public :: &
  pes_x2y1_sysall='x1 x2 y1 x1y1 x2y1'
type (cx_t), public :: &
  pes_x2y1_pc = cx_null
real (kind=dp), allocatable, public :: &
  pes_x2y1_cf(:)

!..procedures
CONTAINS

SUBROUTINE pes_x2y1_read(iun, fn)
integer, intent (in) :: iun
character (len=*), intent (in) :: fn
!-----------------------------------------------------------------------
integer :: nb

open (iun, status='old', file=fn)
read (iun,*) pes_x2y1_pc
read (iun,*) nb
if (nb.ne.pes_x2y1_nb(pes_x2y1_pc%dg)) then
    stop 'pes_x2y1_read: dimension error'
endif
allocate (pes_x2y1_cf(0:nb-1))
if (1.le.nb) then
    read (iun,*) pes_x2y1_cf
endif
close (iun)

return
END SUBROUTINE pes_x2y1_read

FUNCTION pes_x2y1_pot(xn) RESULT (f)
! Potential for generic x2y1
real (kind=dp), intent (in) :: xn(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: r(0:pes_x2y1_nk-1,0:pes_x2y1_nk-1)
integer, parameter :: nki(0:1)=pes_x2y1_nki

call pes_dists(xn, r)
! Dissociates into X1+X1Y1, X2+Y1
f = pes_x1_cf*nki(0)+pes_y1_cf*nki(1)+ &
    cx_f2(nki,r,pes_x2_pc,pes_x2_cf)+ &
    cx_f11(nki,r,pes_x1y1_pc,pes_x1y1_cf)+ &
    cx_f21(nki,r,pes_x2y1_pc,pes_x2y1_cf)

return
END FUNCTION pes_x2y1_pot

FUNCTION pes_x2y1_pot_r(r) RESULT (f)
! Potential for generic x2y1
real (kind=dp), intent(in) :: r(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
integer, parameter :: nki(0:1)=pes_x2y1_nki

! Dissociates into X1+X1Y1, X2+Y1
f = pes_x1_cf*nki(0)+pes_y1_cf*nki(1)+ &
    cx_f2(nki,r,pes_x2_pc,pes_x2_cf)+ &
    cx_f11(nki,r,pes_x1y1_pc,pes_x1y1_cf)+ &
    cx_f21(nki,r,pes_x2y1_pc,pes_x2y1_cf)

return
END FUNCTION pes_x2y1_pot_r

SUBROUTINE pes_x2y1_add(pc, cf)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: nb

if (allocated(pes_x2y1_cf).and.pes_x2y1_pc.eq.pc) then
    pes_x2y1_cf = pes_x2y1_cf+cf
else
    ! Require that pes_x2y1_cf is null
    if (allocated(pes_x2y1_cf)) then
        if (any(pes_x2y1_cf.ne.0.0_dp)) then
            stop 'pes_x2y1_add: mismatch'
        endif
        deallocate (pes_x2y1_cf)
    endif
    nb = pes_x2y1_nb(pc%dg)
    allocate (pes_x2y1_cf(0:nb-1))
    pes_x2y1_pc = pc
    pes_x2y1_cf = cf
endif

return
END SUBROUTINE pes_x2y1_add

SUBROUTINE pes_x2y1_getcf(nki, nkj, pc, cf)
! Callback routine for inv/cx_getcf
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: iun

call pes_getiun(iun)
if (size(nki).ne.size(pes_x2y1_nki)) then
    stop 'pes_x2y1_getcf: bad dimension'
else if (any(nki.ne.pes_x2y1_nki)) then
 stop 'pes_x2y1_getcf: bad nki'
endif
if (all(nkj.eq.(/1,0/))) then
    call pes_x1_add (cf)
    call pes_write0 (iun, 'pcf-x1', pes_x1_cf)
else if (all(nkj.eq.(/0,1/))) then
    call pes_y1_add (cf)
    call pes_write0 (iun, 'pcf-y1', pes_y1_cf)
else if (all(nkj.eq.(/2,0/))) then
    call pes_x2_add (pc, cf)
    call pes_write (iun, 'pcf-x2', pes_x2_pc, pes_x2_cf)
else if (all(nkj.eq.(/1,1/))) then
    call pes_x1y1_add (pc, cf)
    call pes_write (iun, 'pcf-x1y1', pes_x1y1_pc, pes_x1y1_cf)
else if (all(nkj.eq.(/2,1/))) then
    call pes_x2y1_add (pc, cf)
    call pes_write (iun, 'pcf-x2y1', pes_x2y1_pc, pes_x2y1_cf)
else
    stop 'pes_x2y1_getcf: bad nkj'
end if

return
END SUBROUTINE pes_x2y1_getcf

END MODULE pes_x2y1
