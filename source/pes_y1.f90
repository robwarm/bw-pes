MODULE pes_y1
!..use and access
use inv
use pes0
implicit none

private
public :: pes_y1_read, pes_y1_pot, pes_y1_add

!..data
save
integer, parameter, public :: &
  pes_y1_nki(0:0)=(/1/), pes_y1_nk=1, &
  pes_y1_nb=1
character (len=*), parameter, public :: &
  pes_y1_sysall='y1'
real (kind=dp), public :: pes_y1_cf

!..procedures
CONTAINS

SUBROUTINE pes_y1_read(iun, fn)
integer, intent (in) :: iun
character (len=*), intent (in) :: fn
!-----------------------------------------------------------------------

open (iun, status='old', file=fn)
read (iun,*) pes_y1_cf
close (iun)

return
END SUBROUTINE pes_y1_read

FUNCTION pes_y1_pot(xn) RESULT (f)
! Potential for generic Y1
real (kind=dp), intent (in) :: xn(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
! Trivial code, offered for consistency with other pes routines.

f = pes_y1_cf*pes_y1_nki(0)

return
END FUNCTION pes_y1_pot

SUBROUTINE pes_y1_add(cf)
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------

if (size(cf).eq.1) then
    pes_y1_cf = pes_y1_cf+cf(0)
else
    stop 'pes_y1_add: size mismatch'
endif

return
END SUBROUTINE pes_y1_add

END MODULE pes_y1
