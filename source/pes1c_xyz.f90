MODULE pes1c_xyz
!..use and access
use inv
use pes0
use pes1_xyz
implicit none

private
public :: pes1_init_xyz, pes1_init_xyz_consec

CONTAINS

SUBROUTINE pes1_init_xyz(sys)
! Read pcf data files for XYZ systems.
character (len=*), intent (in) :: sys
!------------------------------------------------------------------------
integer :: iun, l0, l1

call pes_getiun(iun)
l0 = len_trim(pes_dir)
l1 = len_trim(sys)
! Read pcf data for the one-body systems
if (cx_substr(sys(1:l1),'x1')) then
    write (*,*) ' reading pcf-x1.dat'
    call pes_x1_read(iun, pes_dir(1:l0)//'pcf-x1.dat')
else
    pes_x1_cf = 0
endif
if (cx_substr(sys(1:l1),'y1')) then
    write (*,*) ' reading pcf-y1.dat'
    call pes_y1_read(iun, pes_dir(1:l0)//'pcf-y1.dat')
else
    pes_y1_cf = 0
endif

! Read pcf data for the multi-body systems.  We've arranged the
! code according to grevlex order on number of X, Y, Z, U.
! Read pcf data for two-body systems.
if (cx_substr(sys(1:l1),'x2')) then
    write (*,*) ' reading pcf-x2.dat'
    call pes_x2_read(iun, pes_dir(1:l0)//'pcf-x2.dat')
else
    pes_x2_pc = cx_null ; allocate (pes_x2_cf(0:-1))
endif
if (cx_substr(sys(1:l1),'x1y1')) then
    write (*,*) ' reading pcf-x1y1.dat'
    call pes_x1y1_read(iun, pes_dir(1:l0)//'pcf-x1y1.dat')
else
    pes_x1y1_pc = cx_null ; allocate (pes_x1y1_cf(0:-1))
endif

! Read pcf data for three-body systems.
if (cx_substr(sys(1:l1),'x2y1')) then
    write (*,*) ' reading pcf-x2y1.dat'
    call pes_x2y1_read(iun, pes_dir(1:l0)//'pcf-x2y1.dat')
else
    pes_x2y1_pc = cx_null ; allocate (pes_x2y1_cf(0:-1))
endif

return
END SUBROUTINE pes1_init_xyz

SUBROUTINE pes1_init_xyz_consec(sys)
! Read pcf data files for XYZ systems.
character (len=*), intent (in) :: sys
!------------------------------------------------------------------------
integer :: iun, l0, l1

call pes_getiun (iun)
l0 = len_trim(pes_dir)
l1 = len_trim(sys)

write (*,*) ' reading consec_pos.dat'
call pes_x2y1consec_read_cp(iun, pes_dir(1:l0)//'consec_pos.dat')
if (cx_substr(sys(1:l1),'x2y1')) then
    write (*,*) ' reading pcf-x2y1consec0.dat'
    call pes_x2y1consec_read(iun, pes_dir(1:l0)//'pcf-x2y1consec0.dat',0)
    write (*,*) ' reading pcf-x2y1consec180.dat'
    call pes_x2y1consec_read(iun, pes_dir(1:l0)//'pcf-x2y1consec180.dat',180)
else
    pes_x2y1consec_pc0 = cx_null_consec ; allocate (pes_x2y1consec_cf0(0:-1))
    pes_x2y1consec_pc180 = cx_null_consec ; allocate (pes_x2y1consec_cf180(0:-1))
endif

return
END SUBROUTINE pes1_init_xyz_consec

END MODULE pes1c_xyz
