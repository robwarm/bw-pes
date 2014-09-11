MODULE pes_x2y1consec
!..use and access
use inv
use pes0
use pes_x2y1
implicit none

private
public :: pes_x2y1consec_read, pes_x2y1consec_read_cp, pes_x2y1consec_pot
public :: pes_x2y1consec_pot_r, pes_x2y1consec_add, pes_x2y1consec_getcf

!..data
save
integer, parameter, public :: &
    pes_x2y1consec_nki(0:1)=(/2,1/), pes_x2y1consec_nk=3, &
    pes_x2y1consec_nb(-1:size(cx_nb21(-1:))-2)=cx_nb21(-1:)
character (len=*), parameter, public :: &
    pes_x2y1consec_sysall='x2y1'
type (cx_t_consec), public :: &
    pes_x2y1consec_pc0 = cx_null_consec, &
    pes_x2y1consec_pc180 = cx_null_consec
type (cp_t_consec), dimension(2), public :: &
    pes_x2y1consec_cp(1:2) = (/ cp_null_consec , cp_null_consec /)
real (kind=dp), allocatable, public :: &
    pes_x2y1consec_cf0(:), pes_x2y1consec_cf180(:)

!..procedures
CONTAINS

SUBROUTINE pes_x2y1consec_read(iun, fn, whichlin)
integer, intent (in) :: iun
character (len=*), intent (in) :: fn
integer, intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: nb

open (iun, status='old', file=fn)

if(whichlin.eq.0) then
    read (iun,*) pes_x2y1consec_pc0
    read (iun,*) nb
    if (nb.ne.pes_x2y1consec_nb(pes_x2y1consec_pc0%dg)) then
        stop 'pes_x2y1consec_read: dimension error'
    endif
    allocate(pes_x2y1consec_cf0(0:nb-1))
    if (1.le.nb) then
        read (iun,*) pes_x2y1consec_cf0
    endif
else
    if(whichlin.eq.180) then
        read (iun,*) pes_x2y1consec_pc180
        read (iun,*) nb
        if (nb.ne.pes_x2y1consec_nb(pes_x2y1consec_pc180%dg)) then
            stop 'pes_x2y1consec_read: dimension error'
        endif
        allocate(pes_x2y1consec_cf180(0:nb-1))
        if (1.le.nb) then
            read (iun,*) pes_x2y1consec_cf180
        endif
    else
        stop 'Illegal whichlin in pes_x2y1consec_read'
    endif
endif

close (iun)

return
END SUBROUTINE pes_x2y1consec_read

SUBROUTINE pes_x2y1consec_read_cp(iun, fn)
integer, intent (in) :: iun
character (len=*), intent (in) :: fn
!-----------------------------------------------------------------------

open (iun, status='old', file=fn)
read (iun,*) pes_x2y1consec_cp(1)
read (iun,*) pes_x2y1consec_cp(2)
close (iun)

return
END SUBROUTINE pes_x2y1consec_read_cp

FUNCTION pes_x2y1consec_pot(xn) RESULT (f)
! Potential for generic x2y1
real (kind=dp), intent (in) :: xn(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: r(0:pes_x2y1consec_nk-1,0:pes_x2y1consec_nk-1)
integer, parameter :: nki(0:1)=pes_x2y1consec_nki

call pes_dists (xn, r)
! Dissociates into X1+X1Y1, X2+Y1

f = pes_x2y1_pot(xn)
f = f + cx_f21_consec(nki,r,pes_x2y1consec_pc0,  pes_x2y1consec_cf0,   pes_x2y1consec_cp, 0) + &
        cx_f21_consec(nki,r,pes_x2y1consec_pc180,pes_x2y1consec_cf180, pes_x2y1consec_cp, 180)

return
END FUNCTION pes_x2y1consec_pot

FUNCTION pes_x2y1consec_pot_r(r) RESULT (f)
! Potential for generic x2y1
real (kind=dp), intent(in) :: r(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------

integer, parameter :: nki(0:1)=pes_x2y1consec_nki

f = pes_x2y1_pot_r(r)
f = f + cx_f21_consec(nki, r, pes_x2y1consec_pc0,   pes_x2y1consec_cf0,   pes_x2y1consec_cp, 0) + &
        cx_f21_consec(nki, r, pes_x2y1consec_pc180, pes_x2y1consec_cf180, pes_x2y1consec_cp, 180)

return
END FUNCTION pes_x2y1consec_pot_r

SUBROUTINE pes_x2y1consec_add(pc, cf, whichlin)
type (cx_t_consec), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
integer, intent (in) :: whichlin
!-----------------------------------------------------------------------
integer :: nb

if(whichlin.eq.0) then
    if (allocated(pes_x2y1consec_cf0).and.pes_x2y1consec_pc0.eq.pc) then
        pes_x2y1consec_cf0 = pes_x2y1consec_cf0+cf
    else
        ! Require that pes_x2y1consec_cf is null
        if (allocated(pes_x2y1consec_cf0)) then
            if (any(dabs(pes_x2y1consec_cf0).ge.1d-12)) then
                stop 'pes_x2y1consec_add: mismatch'
            endif
            deallocate (pes_x2y1consec_cf0)
        endif
        nb = pes_x2y1consec_nb(pc%dg)
        allocate(pes_x2y1consec_cf0(0:nb-1))
        pes_x2y1consec_pc0 = pc
        pes_x2y1consec_cf0 = cf
    endif
else
    if(whichlin.eq.180) then
        if (allocated(pes_x2y1consec_cf180).and.pes_x2y1consec_pc180.eq.pc) then
            pes_x2y1consec_cf180 = pes_x2y1consec_cf180+cf
        else
            ! Require that pes_x2y1consec_cf is null
            if (allocated(pes_x2y1consec_cf180)) then
                if (any(dabs(pes_x2y1consec_cf180).ge.1d-12)) then
                    stop 'pes_x2y1consec_add: mismatch'
                endif
                deallocate (pes_x2y1consec_cf180)
            endif
            nb = pes_x2y1consec_nb(pc%dg)
            allocate (pes_x2y1consec_cf180(0:nb-1))
            pes_x2y1consec_pc180 = pc
            pes_x2y1consec_cf180 = cf
        endif
    else 
        stop 'Illegal whichlin in pes_x2y1consec_add'
    endif
endif

return
END SUBROUTINE pes_x2y1consec_add

SUBROUTINE pes_x2y1consec_getcf(nki, nkj, pc, cf, whichlin)
! Callback routine for inv/cx_getcf
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t_consec), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
integer, intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: iun

call pes_getiun(iun)
if (size(nki).ne.size(pes_x2y1consec_nki)) then
    stop 'pes_x2y1consec_getcf: bad dimension'
else if (any(nki.ne.pes_x2y1_nki)) then
    stop 'pes_x2y1consec_getcf: bad nki'
endif
if (all(nkj.eq.(/2,1/))) then
    call pes_x2y1consec_add (pc, cf,whichlin)
    if(whichlin.eq.0) then
        call pes_write_consec (iun, 'pcf-x2y1consec0', pes_x2y1consec_pc0, pes_x2y1consec_cf0)
    else 
        if(whichlin.eq.180) then
            call pes_write_consec (iun, 'pcf-x2y1consec180', pes_x2y1consec_pc180, pes_x2y1consec_cf180)
        else 
            stop 'Illegal whichlin in pes_x2y1consec_getcf'
        end if
    end if
else
    stop 'pes_x2y1consec_getcf: bad nkj'
end if

return
END SUBROUTINE pes_x2y1consec_getcf

END MODULE pes_x2y1consec
