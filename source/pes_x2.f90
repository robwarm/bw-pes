MODULE pes_x2
!..use and access
use inv
use pes0
use pes_x1
implicit none

private
public :: pes_x2_read, pes_x2_pot, pes_x2_add, pes_x2_getcf

!..data
save
integer, parameter, public :: &
  pes_x2_nki(0:0)=(/2/), pes_x2_nk=2, &
  pes_x2_nb(-1:size(cx_nb2(-1:))-2)=cx_nb2(-1:)
character (len=*), parameter, public :: &
  pes_x2_sysall='x1 x2'
type (cx_t), public :: &
  pes_x2_pc = cx_null
real (kind=dp), allocatable, public :: &
  pes_x2_cf(:)
  
!..procedures
CONTAINS

SUBROUTINE pes_x2_read(iun, fn)
integer, intent (in) :: iun
character (len=*), intent (in) :: fn
!-----------------------------------------------------------------------
integer :: nb

open (iun, status='old', file=fn)
read (iun,*) pes_x2_pc
read (iun,*) nb
if (nb.ne.pes_x2_nb(pes_x2_pc%dg)) then
    stop 'pes_x2_read: dimension error'
endif
allocate (pes_x2_cf(0:nb-1))
if (1.le.nb) then
    read (iun,*) pes_x2_cf
endif
close (iun)

return
END SUBROUTINE pes_x2_read

FUNCTION pes_x2_pot(xn) RESULT (f)
! Potential for generic X2
real (kind=dp), intent (in) :: xn(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
integer, parameter :: nki(0:0)=pes_x2_nki
real (kind=dp) :: r(0:pes_x2_nk-1,0:pes_x2_nk-1)

call pes_dists (xn, r)
f = pes_x1_cf*nki(0)+ &
    cx_f2(nki,r,pes_x2_pc,pes_x2_cf)

return
END FUNCTION pes_x2_pot

SUBROUTINE pes_x2_add(pc, cf)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: nb

if (allocated(pes_x2_cf).and.pes_x2_pc.eq.pc) then
    pes_x2_cf = pes_x2_cf+cf
else
! Require that pes_x2_cf is null
    if (allocated(pes_x2_cf)) then
        if (any(dabs(pes_x2_cf).ge.1d-12)) then
            stop 'pes_x2_add: mismatch'
        endif
        deallocate (pes_x2_cf)
    endif
    nb = pes_x2_nb(pc%dg)
    allocate (pes_x2_cf(0:nb-1))
    pes_x2_pc = pc
    pes_x2_cf = cf
endif

return
END SUBROUTINE pes_x2_add

SUBROUTINE pes_x2_getcf(nki, nkj, pc, cf)
! Callback routine for inv/cx_getcf
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: iun

call pes_getiun (iun)
if (size(nki).ne.size(pes_x2_nki)) then
    stop 'pes_x2_getcf: bad dimension'
else if (any(nki.ne.pes_x2_nki)) then
    stop 'pes_x2_getcf: bad nki'
endif

if (all(nkj.eq.(/1/))) then
    call pes_x1_add (cf)
    call pes_write0 (iun, 'pcf-x1', pes_x1_cf)
else if (all(nkj.eq.(/2/))) then
    call pes_x2_add (pc, cf)
    call pes_write (iun, 'pcf-x2', pes_x2_pc, pes_x2_cf)
else
    stop 'pes_x2_getcf: bad nkj'
end if

return
END SUBROUTINE pes_x2_getcf

END MODULE pes_x2
