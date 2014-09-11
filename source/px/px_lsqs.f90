SUBROUTINE px_lsqs (mat, rhs, rcond, coef)
real (kind=dp), intent (inout) :: mat(0:,0:)
real (kind=dp), intent (in) :: rhs(0:), rcond
real (kind=dp), intent (out) :: coef(0:)
!----------------------------------------------------------------------
integer :: nn, nb, nb0, i, info
real (kind=dp) :: scal(0:size(coef)-1), work(0:4*size(coef)-1), &
  rhs1(0:max(size(coef),size(rhs))-1)
integer :: jpvt(0:size(coef)-1)
nn = size(rhs)
nb = size(coef)
if (size(mat,1).ne.nn.or.size(mat,2).ne.nb) then
 stop 'px_lsqs: bad dimension'
endif
do i = 0, nb-1
 scal(i) = sqrt(sum(mat(0:nn-1,i)**2)/nn)
enddo
write (*,'(1x,a,2x,a,es8.1)') 'Solving the LSQ system: LAPACK dgelsx.', &
  'min(scale) =', minval(scal)
if (any(scal.eq.0)) then
 write (*,'(1x,a)') &
   ' ALERT: some 0 columns; associated coefficients will be 0'
endif
do i = 0, nb-1
 if (scal(i).ne.0) then
  mat(:,i) = mat(:,i)/scal(i)
 endif
enddo
rhs1(0:nn-1) = rhs
jpvt = 0
call dgelsx (nn, nb, 1, mat, nn, rhs1, size(rhs1), jpvt, rcond, &
  nb0, work, info)
! Note: mat has been destroyed
if (info.ne.0) then
 stop 'px_lsqs: LAPACK solver failure'
endif
write (*,'(2x,a,3i6)') &
  'Done dgelsx.  nrhs, ncoef, and rank:', nn, nb, nb0

do i = 0, nb-1
 if (scal(i).eq.0) then
  coef(i) = 0
 else
  coef(i) = rhs1(i)/scal(i)
 endif
enddo
END SUBROUTINE px_lsqs
