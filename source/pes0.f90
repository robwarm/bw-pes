MODULE pes0
!..use and access
use inv
implicit none

private
public :: pes0_init, pes_getiun, pes_nmbrec, pes_dists, &
  pes_write0, pes_write
public :: pes_write_consec

!..data
save
integer, parameter, public :: pes0_dp=selected_real_kind(12,300)
real (kind=dp), parameter, public :: &
  pes_hartree=1.0_dp, pes_joule=pes_hartree/4.359744e-18_dp, &
  pes_wavenmb=pes_hartree/219474.63137_dp, &
  pes_bohr=1.0_dp, pes_angstrom=pes_bohr/0.529177249_dp, &
  pes_elcharge=1.0_dp, pes_coulomb=pes_elcharge/1.60217653e-19_dp
character (len=255), public :: pes_dir='./'
character (len=16), public :: pcf_suffix='.out'

!..procedures
CONTAINS

SUBROUTINE pes0_init(dir, suffix)
! Initialize module data in pes0
character (len=*), intent (in), optional :: dir, suffix
!-----------------------------------------------------------------------
logical :: b0

if (present(dir)) then
    b0 = dir(len_trim(dir):len_trim(dir)).eq.'/'
    if (b0) then
        write (*,*) 'Principal data directory: ', dir(1:len_trim(dir))
        pes_dir = dir
    else
        write (*,*) 'Principal data directory: ', dir(1:len_trim(dir))//'/'
        pes_dir = trim(dir)//'/'
    endif
endif

if (present(suffix)) then
    write (*,*) 'pcf output suffix: ', suffix(1:len_trim(suffix))
    pcf_suffix = suffix
endif

return
END SUBROUTINE pes0_init

SUBROUTINE pes_getiun(iun)
! Obtain a free unit number
integer, intent (out) :: iun
!-----------------------------------------------------------------------
integer :: k
logical :: b

k = 20
inquire (unit=k, opened=b)
do while (b.and.k.lt.100)
    k = k + 1
    inquire (unit=k, opened=b)
enddo

if (.not.b) then
    iun = k
else
    stop 'pes_getiun: no free unit'
endif

return
END SUBROUTINE pes_getiun

FUNCTION pes_nmbrec(iun) RESULT (n)
integer, intent (in) :: iun
integer :: n
!-----------------------------------------------------------------------
integer :: k

k = 0
do while (.true.)
    read (iun,*,end=1)
    k = k + 1
enddo
1 continue
n = k
rewind (iun)

END FUNCTION pes_nmbrec

SUBROUTINE pes_dists(xn, d)
real (kind=dp), intent (in) :: xn(0:,0:)
real (kind=dp), intent (out) :: d(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j, n

if (size(d,1).ne.size(xn,2).or.size(d,2).ne.size(xn,2)) then
    stop 'pes_dists: bad dimensions'
endif

n = size(xn,2)
do j = 0, n-1
    do i = 0, j-1
        d(i,j) = sqrt(sum((xn(:,j)-xn(:,i))**2))
        d(j,i) = d(i,j)
    enddo
    d(j,j) = 0
enddo

return
END SUBROUTINE pes_dists

SUBROUTINE pes_write0(iun, fname, cf0)
integer, intent (in) :: iun
character (len=*), intent (in) :: fname
real (kind=dp), intent (in) :: cf0
!-----------------------------------------------------------------------
integer :: l0, l1

l0 = len_trim(pes_dir)
l1 = len_trim(pcf_suffix)
write (*,*) ' writing ', fname//pcf_suffix(1:l1)
open (iun, status='replace', &
    file=pes_dir(1:l0)//fname//pcf_suffix(1:l1))
write (iun,'(2x,1pg18.12,2x,a)') cf0, 'one-body'
close (iun)

return
END SUBROUTINE pes_write0

SUBROUTINE pes_write(iun, fname, pc, cf)
integer, intent (in) :: iun
character (len=*), intent (in) :: fname
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: l0, l1

l0 = len_trim(pes_dir)
l1 = len_trim(pcf_suffix)
write (*,*) ' writing ', fname//pcf_suffix(1:l1)
open (iun, status='replace', &
    file=pes_dir(1:l0)//fname//pcf_suffix(1:l1))
call cx_write (iun, pc)
write (iun,'(i8,2x,a)') size(cf), 'coefficients'
write (iun,'(4es22.14)') cf
close (iun)

return
END SUBROUTINE pes_write

SUBROUTINE pes_write_consec(iun, fname, pc, cf)
integer, intent (in) :: iun
character (len=*), intent (in) :: fname
type (cx_t_consec), intent (in) :: pc
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: l0, l1

l0 = len_trim(pes_dir)
l1 = len_trim(pcf_suffix)
write (*,*) ' writing ', fname//pcf_suffix(1:l1)
open (iun, status='replace', &
    file=pes_dir(1:l0)//fname//pcf_suffix(1:l1))
call cx_write_consec (iun, pc)
write (iun,'(i8,2x,a)') size(cf), 'coefficients'
write (iun,'(4es22.14)') cf
close (iun)

return
END SUBROUTINE pes_write_consec

END MODULE pes0
