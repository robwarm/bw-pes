SUBROUTINE px_lsq (base, coef)
! Set up and solve least squares system for coef
interface
 subroutine base (xn0, w)
 integer, parameter :: dp=selected_real_kind(12,300)
 real (kind=dp), intent (in) :: xn0(0:,0:)
 real (kind=dp), intent (out) :: w(0:)
 end subroutine base
end interface
real (kind=dp), intent (out) :: coef(0:)
!----------------------------------------------------------------------
integer :: stat0, nd, nk, nb, n, nn0, nn, ng, ng2, ip, i, j, k, &
  i2, j2, ngl, ng2l
real (kind=dp) :: wt, rnd0, rnd1
real (kind=dp), allocatable :: xn1(:,:), vec0(:), vec1(:), vec2(:), &
  vec3(:), rhs(:), mat(:,:)
nd = size(pes2_xn,1)
nk = size(pes2_xn,2)
nb = size(coef)
n = size(pes2_f)
! ng is the number of gradient data that will be used per configuration
if (pes2_havegf) then
 ng = min(px_ng,nk*nd)
else
 ng = 0
endif
! ng2 is the number of hessian data that will be used per configuration
if (pes2_haveg2f) then
 ng2 = min(px_ng2,(nd*nk)*(nd*nk+1)/2)
else
 ng2 = 0
endif
nn0 = 1+ng+ng2
nn = nn0*n
write (*,'(1x,a,i0,a,i0,a,i0,a)') &
  'Forming the LSQ system: ', nn, ' equations, ', &
  nb, ' unknowns; ', nn*nb/1024, ' KW'
allocate (xn1(0:nd-1,0:nk-1), vec0(0:nb-1), vec1(0:nb-1), vec2(0:nb-1), &
  vec3(0:nb-1), rhs(0:nn-1), mat(0:nn-1,0:nb-1), stat=stat0)
if (stat0.ne.0) then
 stop 'px_lsq: failure to allocate workspace'
endif
do ip = 0, n-1
 wt = product(pes2_dwt/(pes2_dwt+pes2_f(ip)-pes2_fmin))
 xn1 = pes2_xn(:,:,ip)
 call base (xn1, vec0)
 mat(ip*nn0,0:nb-1) = wt*vec0
 rhs(ip*nn0) = wt*pes2_df(ip)
 if (0.lt.ng) then
  ngl = 0
  do j = 0, nk-1
   call random_number (rnd0)
   if (rnd0*(nk-j).lt.(ng-ngl+nd-1)/nd) then
    do i = 0, nd-1
     call random_number (rnd1)
     if (rnd1*(nd-i).lt.ng-ngl) then
      xn1 = pes2_xn(:,:,ip)
      xn1(i,j) = xn1(i,j)-pes2_dfdis
      call base (xn1, vec0)
      xn1 = pes2_xn(:,:,ip)
      xn1(i,j) = xn1(i,j)+pes2_dfdis
      call base (xn1, vec1)
      mat(ip*nn0+1+ngl,0:nb-1) = wt*px_wtgf*(vec1-vec0)/(2*pes2_dfdis)
      rhs(ip*nn0+1+ngl) = wt*px_wtgf*pes2_dgf(i,j,ip)
      ngl = ngl+1
     endif
    enddo
   endif
  enddo
  if (ngl.ne.ng) then
   stop 'px_lsq: bad count ngl'
  endif
 endif
 if (0.lt.ng2) then
  ng2l = 0
  k = 0
  do j2 = 0, nd*nk-1
   do i2 = 0, j2
    call random_number (rnd0)
    if (rnd0*((nd*nk)*(nd*nk+1)/2-k).lt.ng2-ng2l) then
     xn1 = pes2_xn(:,:,ip)
     xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)-pes2_df2dis
     xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)-pes2_df2dis
     call base (xn1, vec0)
     xn1 = pes2_xn(:,:,ip)
     xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)+pes2_df2dis
     xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)-pes2_df2dis
     call base (xn1, vec1)
     xn1 = pes2_xn(:,:,ip)
     xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)-pes2_df2dis
     xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)+pes2_df2dis
     call base (xn1, vec2)
     xn1 = pes2_xn(:,:,ip)
     xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)+pes2_df2dis
     xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)+pes2_df2dis
     call base (xn1, vec3)
     mat(ip*nn0+1+ng+ng2l,0:nb-1) = &
       wt*px_wtg2f*(vec0+vec3-vec1-vec2)/(4*pes2_df2dis**2)
     rhs(ip*nn0+1+ng+ng2l) = wt*px_wtg2f*pes2_dg2f(i2,j2,ip)
     ng2l = ng2l+1
    endif
    k = k+1
   enddo
  enddo
  if (ng2l.ne.ng2) then
   stop 'px_lsq: bad count ng2l'
  endif
 endif
enddo
if (1.le.nb) then
 call px_lsqs (mat, rhs, px_rcond, coef)
endif
! Note: mat has been destroyed
!! The following deallocation should not be necessary
deallocate (xn1, vec0, vec1, rhs, mat)
return
END SUBROUTINE px_lsq
