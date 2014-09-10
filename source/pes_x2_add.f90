SUBROUTINE pes_x2_add (pc, cf)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: nb
if (allocated(pes_x2_cf).and.pes_x2_pc.eq.pc) then
 pes_x2_cf = pes_x2_cf+cf
else
! Require that pes_x2_cf is null
 if (allocated(pes_x2_cf)) then
  if (any(pes_x2_cf.ne.0.0_wp)) then
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
