SUBROUTINE pes_x2_getcf (nki, nkj, pc, cf)
! Callback routine for inv/cx_getcf
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: cf(0:)
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
