MODULE pes1c
!..use and access
use inv
use pes1c_xyz
implicit none

private
public :: pes1_init, pes1_init_consec

!..procedures
CONTAINS

SUBROUTINE pes1_init(sys)
! Read pcf data files.
character (len=*), intent (in) :: sys
!------------------------------------------------------------------------

! Announce the intent
write (*,*) 'Expected data: ', sys(1:len_trim(sys))
! Read pcf data for generic xyz systems
call pes1_init_xyz(sys)

return
END SUBROUTINE pes1_init

SUBROUTINE pes1_init_consec(sys)
! Read pcf data files.
character (len=*), intent (in) :: sys
!------------------------------------------------------------------------

! Announce the intent
write (*,*) 'Expected data: ', sys(1:len_trim(sys))
! Read pcf data for generic xyz systems
call pes1_init_xyz_consec(sys)

return
END SUBROUTINE pes1_init_consec

END MODULE pes1c
