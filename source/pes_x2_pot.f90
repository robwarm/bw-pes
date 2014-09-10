FUNCTION pes_x2_pot (xn) RESULT (f)
! Potential for generic X2
real (kind=wp), intent (in) :: xn(0:,0:)
real (kind=wp) :: f
!-----------------------------------------------------------------------
integer, parameter :: nki(0:0)=pes_x2_nki
real (kind=wp) :: r(0:pes_x2_nk-1,0:pes_x2_nk-1)
call pes_dists (xn, r)
f = pes_x1_cf*nki(0)+ &
  cx_f2(nki,r,pes_x2_pc,pes_x2_cf)
return
END FUNCTION pes_x2_pot
