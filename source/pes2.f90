MODULE pes2
!..use and access
use inv
use pes0
implicit none

private
public :: pes2_init, pes2_reinit, pes2_readf, pes2_readf_old, &
  pes2_grad, pes2_grad2, pes2_errf, pes2_errvf, pes2_errf_mod

  !..data
save
integer, parameter, public :: pes2_nd=3
integer, allocatable, public :: pes2_ikord(:), pes2_ikord2(:)
logical, public :: &
  pes2_havegf=.false., pes2_haveg2f=.false., pes2_havedip=.false.
real (kind=dp), public :: &
  pes2_dwt(0:1)=(/0.1_dp,1.0_dp/), &
  pes2_dfdis=1.0e-3_dp, pes2_df2dis=5.0e-3_dp
real (kind=dp), allocatable, public :: &
  pes2_xn(:,:,:), pes2_f(:), pes2_gf(:,:,:), pes2_g2f(:,:,:), &
  pes2_df(:), pes2_dgf(:,:,:), pes2_dg2f(:,:,:), &
  pes2_dip(:,:), pes2_ddip(:,:)
character (len=32), allocatable, public :: &
  pes2_lab(:)
real (kind=dp), public :: &
  pes2_fmin=huge(0.0_dp), pes2_fmax=-huge(0.0_dp)

!..procedures
CONTAINS

SUBROUTINE pes2_init(fname, nk, fun, vfun, iord)
! Initialize module data
integer, intent (in) :: nk
character (len=*), intent (in) :: fname
optional :: vfun
interface
    function fun (xn0) result (f0)
    use pes0, dp=>pes0_dp
    real (kind=dp), intent (in) :: xn0(0:,0:)
    real (kind=dp) :: f0
    end function fun
    function vfun (xn0) result (vf0)
    use pes0, dp=>pes0_dp
    real (kind=dp), intent (in) :: xn0(0:,0:)
    real (kind=dp) :: vf0(0:size(xn0,2)-1)
    end function vfun
end interface
integer, intent (in), optional :: iord(0:nk-1)
!----------------------------------------------------------------------
integer :: iun, n, ip, i, k, iord2(0:pes2_nd*nk-1)
real (kind=dp) :: t0, xn0(0:pes2_nd-1,0:nk-1), xn1(0:pes2_nd-1,0:nk-1), &
  gf0(0:pes2_nd-1,0:nk-1), g2f0(0:pes2_nd*nk-1,0:pes2_nd*nk-1)
character (len=80) :: cht0

call pes_getiun(iun)
open (iun, status='old', file=fname)
n = pes_nmbrec(iun)/(nk+2)
allocate (pes2_ikord(0:nk-1))
allocate (pes2_ikord2(0:pes2_nd*nk-1))
if (present(iord)) then
    pes2_ikord = iord
else
    pes2_ikord = (/(k,k=0,nk-1)/)
endif
do i = 0, pes2_nd*nk-1
    k = i/pes2_nd
    pes2_ikord2(i) = pes2_nd*pes2_ikord(k)+modulo(i,pes2_nd)
enddo
allocate (pes2_xn(0:pes2_nd-1,0:nk-1,0:n-1))
allocate (pes2_lab(0:n-1))
allocate (pes2_f(0:n-1), pes2_df(0:n-1))
if (pes2_havegf) then
    allocate (pes2_gf(0:pes2_nd-1,0:nk-1,0:n-1))
    allocate (pes2_dgf(0:pes2_nd-1,0:nk-1,0:n-1))
endif
if (pes2_haveg2f) then
    allocate (pes2_g2f(0:pes2_nd*nk-1,0:pes2_nd*nk-1,0:n-1))
    allocate (pes2_dg2f(0:pes2_nd*nk-1,0:pes2_nd*nk-1,0:n-1))
endif
if (pes2_havedip) then
    allocate (pes2_dip(0:pes2_nd-1,0:n-1))
    allocate (pes2_ddip(0:pes2_nd-1,0:n-1))
endif
write (*,*) 'Reading ', trim(fname)
cht0 = 'configurations: function values'
if (pes2_havegf) then
    cht0 = trim(cht0)//' and gradients'
endif
if (pes2_haveg2f) then
    cht0 = trim(cht0)//' and hessians'
endif
if (pes2_havedip) then
    cht0 = trim(cht0)//' and dipoles'
endif
write (*,*) 'Expecting', n, trim(cht0)
!! For now, require that pes2_haveg2f implies pes2_havegf
if (pes2_haveg2f.and..not.pes2_havegf) then
    stop 'pes2_init: pes2_haveg2f but not pes2_havegf'
endif
do ip = 0, n-1
    if (pes2_haveg2f) then
        if (pes2_havedip) then
            call pes2_readf(iun, pes2_lab(ip), xn0(:,:), pes2_f(ip), &
                            gf0=gf0(:,:), g2f0=g2f0(:,:), dip=pes2_dip(:,ip))
        else
            call pes2_readf(iun, pes2_lab(ip), xn0(:,:), pes2_f(ip), &
                            gf0=gf0(:,:), g2f0=g2f0(:,:))
        endif
        pes2_xn(:,:,ip) = xn0(:,:)
        pes2_gf(:,:,ip) = gf0(:,:)
        pes2_g2f(:,:,ip) = g2f0(:,:)
    else if (pes2_havegf) then
        if (pes2_havedip) then
            call pes2_readf(iun, pes2_lab(ip), xn0(:,:), pes2_f(ip), &
                            gf0=gf0(:,:), dip=pes2_dip(:,ip))
        else
            call pes2_readf(iun, pes2_lab(ip), xn0(:,:), pes2_f(ip), &
            gf0=gf0(:,:))
        endif
        pes2_xn(:,:,ip) = xn0(:,:)
        pes2_gf(:,:,ip) = gf0(:,:)
    else if (pes2_havedip) then
        call pes2_readf(iun, pes2_lab(ip), xn0(:,:), pes2_f(ip), &
        dip=pes2_dip(:,ip))
        pes2_xn(:,:,ip) = xn0(:,:)
    else
        call pes2_readf(iun, pes2_lab(ip), xn0(:,:), pes2_f(ip))
        pes2_xn(:,:,ip) = xn0(:,:)
    endif
enddo
close (iun)
pes2_fmin = minval(pes2_f)
pes2_fmax = maxval(pes2_f)
write (*,'(2x,a,es14.6,4x,a,es14.6)') ' fmin =', pes2_fmin, ' fmax =', pes2_fmax
do ip = 0, n-1
    xn1 = pes2_xn(:,:,ip)
    pes2_df(ip) = pes2_f(ip)-fun(xn1)
    if (pes2_havegf) then
        pes2_dgf(:,:,ip) = pes2_gf(:,:,ip)- &
        pes2_grad(fun,pes2_dfdis,xn1)
    endif
    if (pes2_haveg2f) then
        pes2_dg2f(:,:,ip) = pes2_g2f(:,:,ip)- &
                            pes2_grad2(fun,pes2_df2dis,xn1)
    endif
    if (pes2_havedip) then
        pes2_ddip(:,ip) = pes2_dip(:,ip)-matmul(xn1,vfun(xn1))
    endif
enddo
write (*,'(2x,a,es14.6,4x,a,es14.6)') 'dfmin =', minval(pes2_df), 'dfmax =', maxval(pes2_df)
if (pes2_havedip) then
    t0 = -huge(t0)
    do ip = 0, n-1
        t0 = max(t0,sqrt(sum(pes2_ddip(:,ip)**2)))
    enddo
    write (*,'(2x,a,es14.6)') 'max norm ddip =', t0
endif

return
END SUBROUTINE pes2_init

SUBROUTINE pes2_reinit(fun, vfun)
! Recompute pes2_df (, pes2_dgf, pes2_dg2f, pes2_ddip)
optional :: vfun
interface
    function fun (xn0) result (f0)
    use pes0, dp=>pes0_dp
    real (kind=dp), intent (in) :: xn0(0:,0:)
    real (kind=dp) :: f0
    end function fun
    function vfun (xn0) result (vf0)
    use pes0, dp=>pes0_dp
    real (kind=dp), intent (in) :: xn0(0:,0:)
    real (kind=dp) :: vf0(0:size(xn0,2)-1)
    end function vfun
end interface
!----------------------------------------------------------------------
integer :: ip

do ip = 0, size(pes2_f)-1
    pes2_df(ip) = pes2_f(ip)-fun(pes2_xn(:,:,ip))
    if (pes2_havegf) then
        pes2_dgf(:,:,ip) = pes2_gf(:,:,ip)- &
                           pes2_grad(fun,pes2_dfdis,pes2_xn(:,:,ip))
    endif
    if (pes2_haveg2f) then
        pes2_dg2f(:,:,ip) = pes2_g2f(:,:,ip)- &
                            pes2_grad2(fun,pes2_df2dis,pes2_xn(:,:,ip))
    endif
    if (pes2_havedip) then
        pes2_ddip(:,ip) = pes2_dip(:,ip)- &
                          matmul(pes2_xn(:,:,ip),vfun(pes2_xn(:,:,ip)))
    endif
enddo

return
END SUBROUTINE pes2_reinit

SUBROUTINE pes2_readf(iun, lab, xn0, f0, gf0, g2f0, dip)
! Read data lab, f0 (, dip), xn0 (, gf0, g2f0)
integer, intent (in) :: iun
character (len=*), intent (out) :: lab
real (kind=dp), intent (out) :: xn0(0:,0:), f0
real (kind=dp), intent (out), optional :: gf0(0:,0:), g2f0(0:,0:), &
  dip(0:)
!-----------------------------------------------------------------------
integer :: nd, nk, k, ik, idum
character (len=8) :: chdum
character (len=1023) :: chrec
real (kind=dp) :: tv0(0:size(xn0,1)-1), tv1(0:size(xn0,1)-1)

nd = size(xn0,1)
nk = size(xn0,2)
if (present(gf0)) then
    if (size(gf0,1).ne.nd.or.size(gf0,2).ne.nk) then
        stop 'pes2_readf: dimension error gf0'
    endif
endif
if (present(g2f0)) then
    if (size(g2f0,1).ne.nd*nk.or.size(g2f0,2).ne.nd*nk) then
        stop 'pes2_readf: dimension error g2f0'
    endif
endif
if (present(dip)) then
    if (size(dip).ne.nd) then
        stop 'pes2_readf: dimension error dip'
    endif
endif
read (iun,*) chrec
read (chrec,*) idum
if (idum.ne.nk) then
    stop 'pes2_readf: bad input nk'
endif
lab = ''
read (chrec,*,end=1) idum, lab
1 continue
! Note the units in our xyz data files: geometry in angstrom,
! energy in hartree, dipole moment in (elcharge)*bohr, energy
! gradient in hartree/bohr, hessian in hartree/bohr^2.
! And note as well: it is the gradient and not the force.
if (present(dip)) then
    read (iun,*) f0, dip
    f0 = f0*pes_hartree
    dip = dip*pes_elcharge*pes_bohr
else
    read (iun,*) f0
    f0 = f0*pes_hartree
endif
do k = 0, nk-1
    if (allocated(pes2_ikord)) then
        ik = pes2_ikord(k)
    else
        ik = k
    endif
    if (present(gf0)) then
        read (iun,*) chdum, tv0(:), tv1(:)
        xn0(:,ik) = tv0(:)*pes_angstrom
        gf0(:,ik) = tv1(:)*pes_hartree/pes_bohr
    else
        read (iun,*) chdum, tv0(:)
        xn0(:,ik) = tv0(:)*pes_angstrom
    endif
enddo
if (present(g2f0)) then
    do k = 0, size(g2f0,2)-1
        if (allocated(pes2_ikord2)) then
            read (iun,*) chdum, g2f0(pes2_ikord2,pes2_ikord2(k))
            g2f0(:,pes2_ikord2(k)) = g2f0(:,pes2_ikord2(k))*pes_hartree/(pes_bohr**2)
        else
            read (iun,*) chdum, g2f0(:,k)
            g2f0(:,k) = g2f0(:,k)*pes_hartree/(pes_bohr**2)
        endif
    enddo
endif

return
END SUBROUTINE pes2_readf

SUBROUTINE pes2_readf_old(iun, xn0, f0, gf0, g2f0)
! Read data xn0, f0 (, gf0)
integer, intent (in) :: iun
real (kind=dp), intent (out) :: xn0(0:,0:), f0
real (kind=dp), intent (out), optional :: gf0(0:,0:), g2f0(0:,0:)
!-----------------------------------------------------------------------
integer :: k

read (iun,*)
do k = 0, size(xn0,2)-1
    read (iun,*) xn0(:,k)
enddo
read (iun,*) f0
if (present(gf0)) then
    do k = 0, size(gf0,2)-1
        read (iun,*) gf0(:,k)
    enddo
endif
if (present(g2f0)) then
    do k = 0, size(g2f0,2)-1
        read (iun,*) g2f0(:,k)
    enddo
endif

return
END SUBROUTINE pes2_readf_old

FUNCTION pes2_grad(fun, dfdis, xn0, acc) RESULT (gf0)
! Compute (grad.fun)(xn0)
interface
    function fun (xn0) result (f0)
    use pes0, dp=>pes0_dp
    real (kind=dp), intent (in) :: xn0(0:,0:)
    real (kind=dp) :: f0
    end function fun
end interface
real (kind=dp), intent (in) :: dfdis, xn0(0:,0:)
integer, intent (in), optional :: acc
real (kind=dp) :: gf0(0:size(xn0,1)-1,0:size(xn0,2)-1)
!-----------------------------------------------------------------------
integer :: iacc, i, j
real (kind=dp) :: fa, fb, fc, fd, xn1(0:size(xn0,1)-1,0:size(xn0,2)-1)

if (present(acc)) then
    iacc = acc
else
    iacc = 0
endif
do j = 0, size(xn0,2)-1
    do i = 0, size(xn0,1)-1
        xn1 = xn0
        xn1(i,j) = xn0(i,j)-dfdis
        fa = fun(xn1)
        xn1(i,j) = xn0(i,j)+dfdis
        fb = fun(xn1)
        if (iacc.eq.0) then
            gf0(i,j) = (fb-fa)/(2*dfdis)
        else
            xn1(i,j) = xn0(i,j)-2*dfdis
            fc = fun(xn1)
            xn1(i,j) = xn0(i,j)+2*dfdis
            fd = fun(xn1)
            gf0(i,j) = (8*(fb-fa)-(fd-fc))/(12*dfdis)
        endif
    enddo
enddo

return
END FUNCTION pes2_grad

FUNCTION pes2_grad2(fun, df2dis, xn0) RESULT (g2f0)
! Compute (grad.grad.fun)(xn0)
interface
    function fun (xn0) result (f0)
    use pes0, dp=>pes0_dp
    real (kind=dp), intent (in) :: xn0(0:,0:)
    real (kind=dp) :: f0
    end function fun
end interface
real (kind=dp), intent (in) :: df2dis, xn0(0:,0:)
real (kind=dp) :: g2f0(0:size(xn0)-1,0:size(xn0)-1)
!-----------------------------------------------------------------------
integer :: nd, i2, j2
real (kind=dp) :: fa, fb, fc, fd, xn1(0:size(xn0,1)-1,0:size(xn0,2)-1)

nd = size(xn0,1)
do j2 = 0, size(xn0)-1
    do i2 = 0, j2-1
        xn1 = xn0
        xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)-df2dis
        xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)-df2dis
        fa = fun(xn1)
        xn1 = xn0
        xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)+df2dis
        xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)-df2dis
        fb = fun(xn1)
        xn1 = xn0
        xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)-df2dis
        xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)+df2dis
        fc = fun(xn1)
        xn1 = xn0
        xn1(modulo(i2,nd),i2/nd) = xn1(modulo(i2,nd),i2/nd)+df2dis
        xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)+df2dis
        fd = fun(xn1)
        g2f0(i2,j2) = (fa+fd-fb-fc)/(4*df2dis**2)
        g2f0(j2,i2) = g2f0(i2,j2)
    enddo
    xn1 = xn0
    xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)-df2dis
    fa = fun(xn1)
    xn1 = xn0
    fb = fun(xn1)
    xn1 = xn0
    xn1(modulo(j2,nd),j2/nd) = xn1(modulo(j2,nd),j2/nd)+df2dis
    fc = fun(xn1)
    g2f0(j2,j2) = (fa+fc-2*fb)/(df2dis**2)
    enddo

return
END FUNCTION pes2_grad2

SUBROUTINE pes2_errf(nk, lev)
! Print an error report
integer, intent (in) :: nk, lev
!-----------------------------------------------------------------------
integer :: n, n1, n2, n5, ip, k
real (kind=dp) :: wt, t0, tw, t1, t2, t5, t1d, t2d, t5d, &
  err, errw, err1, err2, err5, errd1, errd2, errd5
real (kind=dp) :: gf0(0:pes2_nd-1,0:nk-1)

t0 = 0 ; tw = 0
t1 = 0 ; t2 = 0 ; t5 = 0
t1d = 0 ; t2d = 0 ; t5d = 0
n1 = 0 ; n2 = 0 ; n5 = 0
n = size(pes2_f)
do ip = 0, n-1
    wt = product(pes2_dwt/(pes2_dwt+pes2_f(ip)-pes2_fmin))
    t0 = t0+pes2_df(ip)**2
    tw = tw+(wt*pes2_df(ip))**2
    if (pes2_f(ip)-pes2_fmin.lt.pes2_dwt(0)) then
        t1 = t1+pes2_df(ip)**2
        if (pes2_havedip) then
            t1d = t1d+sum(pes2_ddip(:,ip)**2)
        endif
        n1 = n1+1
    else if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0)) then
        t2 = t2+pes2_df(ip)**2
        if (pes2_havedip) then
            t2d = t2d+sum(pes2_ddip(:,ip)**2)
        endif
        n2 = n2+1
    else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0)) then
        t5 = t5+pes2_df(ip)**2
        if (pes2_havedip) then
            t5d = t5d+sum(pes2_ddip(:,ip)**2)
        endif
        n5 = n5+1
    endif
enddo
err = sqrt(t0/n)
errw = sqrt(tw/n)
err1 = sqrt(t1/max(1,n1))
err2 = sqrt(t2/max(1,n2))
err5 = sqrt(t5/max(1,n5))
if (pes2_havedip) then
    errd1 = sqrt(t1d/max(1,n1))
    errd2 = sqrt(t2d/max(1,n2))
    errd5 = sqrt(t5d/max(1,n5))
endif
write (*,'(2x,i0,1x,a,1p,(4g9.2))') &
        n, 'configurations, weighting parameter dwt(0:) =', pes2_dwt
write (*,'(2x,es9.2,2x,a,4x,a,es9.2)') &
        err, 'rms error over all configs;', &
        'weighted rms =', errw
write (*,'(2x,a)') &
        'number and rms error for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
write (*,'(2x,3(2x,i7,1x,es10.2))') &
        n1, err1, n2, err2, n5, err5
if (pes2_havedip) then
    write (*,'(2x,a)') &
        'dipole errors for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
    write (*,'(2x,3(2x,i7,1x,es10.2))') &
        n1, errd1, n2, errd2, n5, errd5
endif
if (1.le.lev) then
    ! look for outliers on pes2_f
    write (*,'(a)') ' Looking for outliers on f.'
    do ip = 0, n-1
        if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0).and. &
                    5*max(err1,err2).lt.abs(pes2_df(ip))) then
            write (*,'(2x,i8,3f12.4,f10.2,1x,a)') &
                    ip, pes2_f(ip), pes2_f(ip)-pes2_fmin, pes2_df(ip), &
                    pes2_df(ip)/max(err1,err2), '(*err12)'
        else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0).and. &
                    5*max(err1,err2,err5).lt.abs(pes2_df(ip))) then
            write (*,'(2x,i8,3f12.4,f10.2,1x,a)') &
                    ip, pes2_f(ip), pes2_f(ip)-pes2_fmin, pes2_df(ip), &
                    pes2_df(ip)/max(err1,err2,err5), &
                    '(*err125)'
        endif
        if (1.le.ip) then
        !!    if (pes2_f(ip-1).eq.pes2_f(ip)) then
        !!     write (*,'(a,2x,i8)') ' Alert, equal f:', ip
        !!    endif
        endif
    enddo
endif
if (2.le.lev.and.pes2_havegf) then
    ! look for extreme outliers on pes2_gf
    write (*,'(a)') ' Looking for extreme outliers on grad.f.'
    do ip = 0, n-1
        gf0 = pes2_gf(:,:,ip)-pes2_dgf(:,:,ip)
        t0 = sqrt(sum(gf0(:,:)**2))
        t1 = sqrt(sum(pes2_dgf(:,:,ip)**2))
        if (t0.lt.t1) then
            write (*,'(a,i6)') ' ip =', ip
            !! reorder the nuclei?
            do k = 0, nk-1
                write (*,'(a6,i4,4f16.8)') 'gf', k, pes2_gf(:,k,ip)
                write (*,'(a6,i4,4f16.8)') 'gf0', k, gf0(:,k)
                write (*,'(a6,i4,4f16.8)') 'dgf', k, pes2_dgf(:,k,ip)
            enddo
        endif
    enddo
endif
if (1.le.lev.and.pes2_havedip) then
    ! look for outliers on pes2_dip
    write (*,'(a)') ' Looking for outliers on dip.'
    do ip = 0, n-1
        t0 = sqrt(sum(pes2_ddip(:,ip)**2))
        if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0).and. &
                        5*max(errd1,errd2).lt.t0) then
            write (*,'(2x,i8,4f12.4,f10.2,1x,a)') &
                        ip, pes2_f(ip), pes2_dip(:,ip), &
                        t0/max(errd1,errd2), '(*errd12)'
        else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0).and. &
                        5*max(errd1,errd2,errd5).lt.t0) then
            write (*,'(2x,i8,4f12.4,f10.2,1x,a)') &
                        ip, pes2_f(ip), pes2_dip(:,ip), &
                        t0/max(errd1,errd2,errd5), &
                        '(*errd125)'
        endif
    enddo
endif

return
END SUBROUTINE pes2_errf

SUBROUTINE pes2_errvf(nk, lev)
! Print an error report
integer, intent (in) :: nk, lev
!-----------------------------------------------------------------------
integer :: n, n1, n2, n5, ip, k
real (kind=dp) :: wt, t0, t1d, t2d, t5d, &
  err, err1, err2, err5, errd1, errd2, errd5
real (kind=dp) :: gf0(0:pes2_nd-1,0:nk-1)

t1d = 0 ; t2d = 0 ; t5d = 0
n1 = 0 ; n2 = 0 ; n5 = 0
n = size(pes2_f)
if (pes2_havedip) then
    do ip = 0, n-1
        wt = product(pes2_dwt/(pes2_dwt+pes2_f(ip)-pes2_fmin))
        if (pes2_f(ip)-pes2_fmin.lt.pes2_dwt(0)) then
   t1d = t1d+sum(pes2_ddip(:,ip)**2)
   n1 = n1+1
  else if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0)) then
   t2d = t2d+sum(pes2_ddip(:,ip)**2)
   n2 = n2+1
  else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0)) then
   t5d = t5d+sum(pes2_ddip(:,ip)**2)
   n5 = n5+1
  endif
 enddo
 errd1 = sqrt(t1d/max(1,n1))
 errd2 = sqrt(t2d/max(1,n2))
 errd5 = sqrt(t5d/max(1,n5))
 write (*,'(2x,i0,1x,a,1p,(4g9.2))') &
   n, 'configurations, weighting parameter dwt(0:) =', pes2_dwt
 write (*,'(2x,a)') &
   'number and rms error for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
 write (*,'(2x,a)') &
   'dipole errors for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
 write (*,'(2x,3(2x,i7,1x,es10.2))') &
   n1, errd1, n2, errd2, n5, errd5
 if (1.le.lev) then
! look for outliers on pes2_dip
  write (*,'(a)') ' Looking for outliers on dip.'
  do ip = 0, n-1
   t0 = sqrt(sum(pes2_ddip(:,ip)**2))
   if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0).and. &
     5*max(errd1,errd2).lt.t0) then
    write (*,'(2x,i8,4f12.4,f10.2,1x,a)') &
      ip, pes2_f(ip), pes2_dip(:,ip), &
      t0/max(errd1,errd2), '(*errd12)'
   else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0).and. &
     5*max(errd1,errd2,errd5).lt.t0) then
    write (*,'(2x,i8,4f12.4,f10.2,1x,a)') &
      ip, pes2_f(ip), pes2_dip(:,ip), &
      t0/max(errd1,errd2,errd5), &
      '(*errd125)'
   endif
  enddo
 endif
endif
return
END SUBROUTINE pes2_errvf

SUBROUTINE pes2_errf_mod (nk, lev)
! Print an error report
integer, intent (in) :: nk, lev
!-----------------------------------------------------------------------
real (kind=dp), parameter :: pi=3.14159265358979323846_dp
integer :: n, n1, n2, n5, ip, k,nmod,nmod2,nmod3
real (kind=dp) :: wt, t0, tw, t1, t2, t5, t1d, t2d, t5d, &
  err, errw, err1, err2, err5, errd1, errd2, errd5, errmod,errmod2,errmod3,tmod,tmod2,tmod3
real (kind=dp) :: gf0(0:pes2_nd-1,0:nk-1)

real (kind=dp) :: s,a,b,c,d,alpha, xn1(0:2,0:nk-1),r(0:nk-1,0:nk-1)
t0 = 0 ; tw = 0
t1 = 0 ; t2 = 0 ; t5 = 0
t1d = 0 ; t2d = 0 ; t5d = 0
n1 = 0 ; n2 = 0 ; n5 = 0
tmod=0.0_dp

n = size(pes2_f)
do ip = 0, n-1
 xn1 = pes2_xn(:,:,ip)
call pes_dists(xn1,r)
a = r(1,0) !r_HH
b = r(2,0) !r_HC
c = r(2,1) !R_HC
s=(a+b+c)/2._dp
d=dsqrt( (s-a)*(s-b)*(s-c)/s )
alpha= 2._dp*atan( d/(s-a))
if(isnan(alpha)) alpha=pi

 wt = product(pes2_dwt/(pes2_dwt+pes2_f(ip)-pes2_fmin))
 t0 = t0+pes2_df(ip)**2
 tw = tw+(wt*pes2_df(ip))**2
 if (pes2_f(ip)-pes2_fmin.lt.pes2_dwt(0)) then
  t1 = t1+pes2_df(ip)**2
  if (pes2_havedip) then
   t1d = t1d+sum(pes2_ddip(:,ip)**2)
  endif
  n1 = n1+1
 else if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0)) then
  t2 = t2+pes2_df(ip)**2
  if (pes2_havedip) then
   t2d = t2d+sum(pes2_ddip(:,ip)**2)
  endif
  n2 = n2+1
 else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0)) then
  t5 = t5+pes2_df(ip)**2
  if (pes2_havedip) then
   t5d = t5d+sum(pes2_ddip(:,ip)**2)
  endif
  n5 = n5+1
 endif

 if ( (pes2_f(ip)-pes2_fmin.lt.0.185_dp).and.(pes2_f(ip)-pes2_fmin.gt.0.15_dp)) then
  tmod = tmod+pes2_df(ip)**2
  !if (pes2_havedip) then
  ! t1d = t1d+sum(pes2_ddip(:,ip)**2)
  !endif
  nmod=nmod+1
 end if
 if ( (pes2_f(ip)-pes2_fmin.lt.0.185_dp).and.(pes2_f(ip)-pes2_fmin.gt.0.15_dp).and.(alpha.le.(pi/4._dp)).and.(max(b,c).lt.16.0_dp ) ) then
  tmod2 = tmod2+pes2_df(ip)**2
  !if (pes2_havedip) then
  ! t1d = t1d+sum(pes2_ddip(:,ip)**2)
  !endif
  nmod2=nmod2+1
 end if
 if ( (pes2_f(ip)-pes2_fmin.lt.0.185_dp).and.(pes2_f(ip)-pes2_fmin.gt.0.15_dp).and.(alpha.ge.(3._dp*pi/4._dp)).and.(max(b,c).lt.16.0_dp ) ) then
  tmod3 = tmod3+pes2_df(ip)**2
  !if (pes2_havedip) then
  ! t1d = t1d+sum(pes2_ddip(:,ip)**2)
  !endif
  nmod3=nmod3+1
 end if



enddo
err = sqrt(t0/n)
errw = sqrt(tw/n)
err1 = sqrt(t1/max(1,n1))
err2 = sqrt(t2/max(1,n2))
err5 = sqrt(t5/max(1,n5))

errmod=sqrt(tmod/max(1,nmod))
errmod2=sqrt(tmod2/max(1,nmod2))
errmod3=sqrt(tmod3/max(1,nmod3))
if (pes2_havedip) then
 errd1 = sqrt(t1d/max(1,n1))
 errd2 = sqrt(t2d/max(1,n2))
 errd5 = sqrt(t5d/max(1,n5))
endif
write (*,'(2x,i0,1x,a,1p,(4g9.2))') &
  n, 'configurations, weighting parameter dwt(0:) =', pes2_dwt
write (*,'(2x,es9.2,2x,a,4x,a,es9.2)') &
  err, 'rms error over all configs;', &
  'weighted rms =', errw
write (*,'(2x,a)') &
  'number and rms error for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
write (*,'(2x,3(2x,i7,1x,es10.2))') &
  n1, err1, n2, err2, n5, err5

write(555,'(8es16.5)') err1,err2,err5,err, errw, errmod,errmod2,errmod3

if (pes2_havedip) then
 write (*,'(2x,a)') &
   'dipole errors for (f0(ip)-fmin)/dwt(0) in [0,1), [1,2), [2,5):'
 write (*,'(2x,3(2x,i7,1x,es10.2))') &
   n1, errd1, n2, errd2, n5, errd5
endif
if (1.le.lev) then
! look for outliers on pes2_f
 write (*,'(a)') ' Looking for outliers on f.'
 do ip = 0, n-1
  if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0).and. &
    5*max(err1,err2).lt.abs(pes2_df(ip))) then
   write (*,'(2x,i8,3f12.4,f10.2,1x,a)') &
     ip, pes2_f(ip), pes2_f(ip)-pes2_fmin, pes2_df(ip), &
     pes2_df(ip)/max(err1,err2), '(*err12)'
  else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0).and. &
    5*max(err1,err2,err5).lt.abs(pes2_df(ip))) then
   write (*,'(2x,i8,3f12.4,f10.2,1x,a)') &
     ip, pes2_f(ip), pes2_f(ip)-pes2_fmin, pes2_df(ip), &
     pes2_df(ip)/max(err1,err2,err5), &
     '(*err125)'
  endif
  if (1.le.ip) then
!!    if (pes2_f(ip-1).eq.pes2_f(ip)) then
!!     write (*,'(a,2x,i8)') ' Alert, equal f:', ip
!!    endif
  endif
 enddo
endif
if (2.le.lev.and.pes2_havegf) then
! look for extreme outliers on pes2_gf
 write (*,'(a)') ' Looking for extreme outliers on grad.f.'
 do ip = 0, n-1
  gf0 = pes2_gf(:,:,ip)-pes2_dgf(:,:,ip)
  t0 = sqrt(sum(gf0(:,:)**2))
  t1 = sqrt(sum(pes2_dgf(:,:,ip)**2))
  if (t0.lt.t1) then
   write (*,'(a,i6)') ' ip =', ip
!! reorder the nuclei?
   do k = 0, nk-1
    write (*,'(a6,i4,4f16.8)') 'gf', k, pes2_gf(:,k,ip)
    write (*,'(a6,i4,4f16.8)') 'gf0', k, gf0(:,k)
    write (*,'(a6,i4,4f16.8)') 'dgf', k, pes2_dgf(:,k,ip)
   enddo
  endif
 enddo
endif
if (1.le.lev.and.pes2_havedip) then
! look for outliers on pes2_dip
 write (*,'(a)') ' Looking for outliers on dip.'
 do ip = 0, n-1
  t0 = sqrt(sum(pes2_ddip(:,ip)**2))
  if (pes2_f(ip)-pes2_fmin.lt.2*pes2_dwt(0).and. &
    5*max(errd1,errd2).lt.t0) then
   write (*,'(2x,i8,4f12.4,f10.2,1x,a)') &
     ip, pes2_f(ip), pes2_dip(:,ip), &
     t0/max(errd1,errd2), '(*errd12)'
  else if (pes2_f(ip)-pes2_fmin.lt.5*pes2_dwt(0).and. &
    5*max(errd1,errd2,errd5).lt.t0) then
   write (*,'(2x,i8,4f12.4,f10.2,1x,a)') &
     ip, pes2_f(ip), pes2_dip(:,ip), &
     t0/max(errd1,errd2,errd5), &
     '(*errd125)'
  endif
 enddo
endif
return
END SUBROUTINE pes2_errf_mod

END MODULE pes2
