MODULE pes_x1
!..use and access
use inv
use pes0
implicit none

private
public :: pes_x1_read, pes_x1_pot, pes_x1_add

!..data
save
integer, parameter, public :: &
  pes_x1_nki(0:0)=(/1/), pes_x1_nk=1, &
  pes_x1_nb=1
character (len=*), parameter, public :: &
  pes_x1_sysall='x1'
real (kind=dp), public :: pes_x1_cf

!..procedures
CONTAINS

SUBROUTINE pes_x1_read(iun, fn)
integer, intent (in) :: iun
character (len=*), intent (in) :: fn
!-----------------------------------------------------------------------

open (iun, status='old', file=fn)
read (iun,*) pes_x1_cf
close (iun)

return
END SUBROUTINE pes_x1_read

FUNCTION pes_x1_pot(xn) RESULT (f)
! Potential for generic X1
real (kind=dp), intent (in) :: xn(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
! Trivial code, offered for consistency with other pes routines.

f = pes_x1_cf*pes_x1_nki(0)

return
END FUNCTION pes_x1_pot

SUBROUTINE pes_x1_add(cf)
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------

if (size(cf).eq.1) then
    pes_x1_cf = pes_x1_cf+cf(0)
else
    stop 'pes_x1_add: size mismatch'
endif

return
END SUBROUTINE pes_x1_add

END MODULE pes_x1
