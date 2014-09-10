MODULE inv_any
!..use and access
implicit none

private
public :: any_base, any_mkprib

!..procedures
CONTAINS

SUBROUTINE any_base(mxd, npr, nsc, nb, u, v, w)
integer, intent (in) :: mxd, npr(0:), nsc(0:), nb(0:)
real (kind=dp), intent (in) :: u(0:), v(0:)
real (kind=dp), intent (out) :: w(0:nb(mxd)-1)
!-----------------------------------------------------------------------
integer :: l(0:mxd,0:npr(mxd)), ind, i, k, d, inc

if (size(u).lt.npr(mxd).or.size(v).lt.nsc(mxd)) then
    stop 'any_base: bad dimensions'
endif

l(0,0:npr(mxd)) = 0
if (1.le.size(w)) then
    w(0) = v(0)
endif
ind = 1
do d = 1, mxd
    do k = 1, d
        do i = npr(k-1), npr(k)-1
            l(d,i) = ind
            inc = l(d+1-k,0) - l(d-k,i)
            w(ind:ind+inc-1) = u(i) * w(l(d-k,i):l(d-k,i)+inc-1)
            ind = ind + inc
        enddo
    enddo
    l(d,npr(d):npr(mxd)) = ind
    w(ind:ind+nsc(d)-nsc(d-1)-1) = v(nsc(d-1):nsc(d)-1)
    ind = ind+nsc(d)-nsc(d-1)
enddo

END SUBROUTINE any_base

SUBROUTINE any_mkprib(name, mxd, dnpr, dnpb)
character (len=*), intent (in) :: name
integer, intent (in) :: mxd, dnpr(0:), dnpb(0:)
!-----------------------------------------------------------------------
integer, allocatable :: l(:,:)
integer :: ind, npr, i, k, d, inc
character (len=*), parameter :: fmt0='(a)', &
  fmtn='(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)'

allocate (l(0:mxd,0:sum(dnpr(0:mxd))))
write (*,fmt0) &
  'SUBROUTINE '//name//'_prib (mxd, u, w)', &
  'integer, intent (in) :: mxd', &
  'real (kind=dp), intent (in) :: u(0:)', &
  'real (kind=dp), intent (out) :: w(0:'//name//'_npb(mxd)-1)', &
  '!---------------------------------------'// &
    '--------------------------------', &
  'if (size(u).ne.nr) then', &
  " stop '"//name//"_prib: bad size u'", &
  '!! else if (size(w).ne.'//name//'_npb(mxd)) then', &
  "!!  stop '"//name//"_prib: bad size w'", &
  'endif'
write (*,fmt0) &
  '! The following code was obtained using these parameters'
write (*,'(a,10(1x,i0):(" &"/"!  ",10(1x,i0)))') &
  '! prims, dnpr(0:*):', dnpr(0:mxd)
write (*,'(a,10(1x,i0):(" &"/"!  ",10(1x,i0)))') &
  '! prib,  dnpb(0:*):', dnpb(0:mxd)
l(0,:) = 0
write (*,fmt0) &
  '! constant term', &
  'w(0) = 1'
  
ind = 1
do d = 1, mxd
    write (*,fmtn) '! terms of degree ', d
    write (*,fmtn) 'if (', d, '.le.mxd) then'
    npr = dnpr(0)
    do k = 1, d
        do i = npr, npr+dnpr(k)-1
            l(d,i) = ind
            inc = l(d+1-k,0) - l(d-k,i)
            if (inc.eq.1) then
                write (*,fmtn) ' w(', ind, ') = u(', i, ')*w(', l(d-k,i), ')'
            else if (2.le.inc) then
                write (*,fmtn) ' w(', ind, ':', ind+inc-1, ') = u(', i, ')*w(', &
                                l(d-k,i), ':', l(d-k,i)+inc-1, ')'
            endif
            ind = ind + inc
        enddo
        npr = npr + dnpr(k)
    enddo
    l(d,npr:sum(dnpr(0:mxd))) = ind
    write (*,fmt0) 'endif'
enddo
write (*,fmt0) 'END SUBROUTINE '//name//'_prib'

END SUBROUTINE any_mkprib

END MODULE inv_any
