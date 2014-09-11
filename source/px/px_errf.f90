SUBROUTINE px_errf (base, coef)
! Print an error report
interface
 subroutine base (xn0, w)
 integer, parameter :: dp=selected_real_kind(12,300)
 real (kind=dp), intent (in) :: xn0(0:,0:)
 real (kind=dp), intent (out) :: w(0:)
 end subroutine base
end interface
real (kind=dp), intent (in) :: coef(0:)
!-----------------------------------------------------------------------
integer :: n, n1, n2, n5, ip
real (kind=dp) :: wt, df0, t0, tw, t1, t2, t5, err, errw, &
  err1, err2, err5
real (kind=dp) :: vec(0:size(coef)-1)
t0 = 0 ; tw = 0
t1 = 0 ; t2 = 0 ; t5 = 0
n1 = 0 ; n2 = 0 ; n5 = 0
n = size(pes2_f)
do ip = 0, n-1
 wt = product(pes2_dwt/(pes2_dwt+pes2_f(ip)-pes2_fmin))
 call base (pes2_xn(:,:,ip), vec)
 df0 = pes2_df(ip)-dot_product(coef,vec)
 t0 = t0+df0**2
 tw = tw+(wt*df0)**2
 if (pes2_f(ip)-pes2_fmin.lt.pes2_dwt(0)) then
  t1 = t1+df0**2
  n1 = n1+1
 else if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0)) then
  t2 = t2+df0**2
  n2 = n2+1
 else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0)) then
  t5 = t5+df0**2
  n5 = n5+1
 endif
enddo
err = sqrt(t0/n)
errw = sqrt(tw/n)
err1 = sqrt(t1/max(1,n1))
err2 = sqrt(t2/max(1,n2))
err5 = sqrt(t5/max(1,n5))
write (*,'(2x,i0,1x,a,1p,(4g9.2))') &
  n, 'configurations, weighting parameter dwt(0:) =', pes2_dwt
write (*,'(2x,es9.2,2x,a,4x,a,es9.2)') &
  err, 'rms error over all configs;', &
  'weighted rms =', errw
write (*,'(2x,a)') &
  'number and rms error for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
write (*,'(2x,3(2x,i7,1x,es10.2))') &
  n1, err1, n2, err2, n5, err5
return
END SUBROUTINE px_errf
