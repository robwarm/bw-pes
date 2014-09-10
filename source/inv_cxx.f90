MODULE inv_cxx
!..use and access
use inv_core
implicit none

private
public :: cx_cut, cx_var, cx_write, operator(.eq.), cx_sort, cx_dim, &
  cx_substr, cx_charid

public :: cx_cut_consec, cx_var_consec, cx_write_consec
!..types
type, public :: cx_t
 integer :: dg, kx, lx
 real (kind=dp) :: ax, bx
end type

!dg-degree of polynomial
!ax,kx,lx- as usual for cx_cut_consec
!r0x, bx,cx for cx_var_consec
type, public :: cx_t_consec
 integer :: dg
 real (kind=dp) :: ax,kx,lx, r0x,bx,cx
end type

! r0x = (a * alpha + b) * r + c * alpha + d
type, public :: cp_t_consec
 real (kind=dp) :: a, b, c, d
end type

!..data
type (cx_t), parameter, public :: cx_null=cx_t(-1,0,0,0.0_dp,0.0_dp)
type (cx_t_consec), parameter, public :: cx_null_consec=cx_t_consec(-1,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp)
type (cp_t_consec), parameter, public :: cp_null_consec=cp_t_consec(0.0_dp,0.0_dp,0.0_dp,0.0_dp)

!..procedures
interface operator (.eq.)
    module procedure cx_equals
end interface

interface operator (.eq.)
    module procedure cx_equals_consec
end interface

CONTAINS

FUNCTION cx_cut(pc, r) result (f)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
integer :: i, j
real (kind=dp) :: d2, d

if (size(r,1).ne.size(r,2).or.size(r,1).lt.2) then
    stop 'cx_cut: bad dimensions'
endif

d2 = 0
do j = 1, size(r,2)-1
    do i = 0, j-1
        d2 = d2 + r(i,j)**2 + r(j,i)**2
    enddo
enddo

if (d2.eq.0) then
    stop 'cx_cut: zero distance'
endif

d = sqrt(d2/(size(r,2)*(size(r,2)-1)))
if (d/pc%ax.lt.1) then
    f = (1-d/pc%ax)**pc%kx/d**pc%lx
else
    f = 0
endif

return
END FUNCTION cx_cut

FUNCTION cx_cut_consec(pc, r,whichlin) result (f)
type (cx_t_consec), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp) :: f
integer,intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: i, j
real (kind=dp) :: d2, d
real (kind=dp) :: s,a,c,b, alpha
real (kind=dp), parameter :: pi=3.14159265358979323846_dp

if (size(r,1).ne.size(r,2).or.size(r,1).lt.2) then
    stop 'cx_cut_consec: bad dimensions'
endif

d2 = 0
do j = 1, size(r,2)-1
    do i = 0, j-1
        d2 = d2+r(i,j)**2+r(j,i)**2
    enddo
enddo

if (d2.eq.0) then
 stop 'cx_cut_consec: zero distance'
endif

d = sqrt(d2/(size(r,2)*(size(r,2)-1)))

!for concical intersection: function should be localized around intersection
!                               --> e.g. exp(-alpha) * exp(- (r-r_0)^2) 

!irregular triangle:
!s=(a+b+c)/2
!r=sqrt( (s-a)(s-b)(s-c)/s )
!alpha= 2 atan( r/s-a)

a = r(1,0) !r_HH
b = r(2,0) !r_HC
c = r(2,1) !R_HC

s=(a+b+c)/2._dp
d=dsqrt( (s-a)*(s-b)*(s-c)/s )
alpha= 2._dp*atan( d/(s-a))
if(isnan(alpha)) alpha=pi

if (whichlin.eq.0) then
!ok, the C-H-H case: the angle H-C-H is zero :
    f = exp(-pc%kx*alpha)*max(exp(-pc%lx*dabs(b-pc%ax)),exp(-pc%lx*dabs( c-pc%ax)))
else if(whichlin.eq.180) then
    f = exp(-pc%kx*(pi-alpha))*max(exp(-pc%lx*dabs(b-pc%ax)),exp(-pc%lx*dabs( c-pc%ax)))
else
    stop 'Illegal whichlin in cx_cut_consec'
end if

return
END FUNCTION cx_cut_consec

SUBROUTINE cx_var(pc, r, y)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: y(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j

if (size(r,1).ne.size(r,2).or.size(y,1).ne.size(r,1).or. &
    size(y,2).ne.size(r,2)) then
    stop 'cx_var: bad dimensions'
endif

do j = 0, size(r,2)-1
    do i = 0, size(r,1)-1
        if (i.ne.j) then
            if (pc%bx.eq.0.0_dp) then
                y(i,j) = r(i,j)
            else if (pc%bx.gt.0) then
            !  exponential transformation
                y(i,j) = exp(-r(i,j)/pc%bx)
            else if (pc%bx.lt.0) then
            ! linear at small r, tends to abs(pc%bx) at large r
                y(i,j) = r(i,j)/sqrt(1+(r(i,j)/pc%bx)**2)
            endif
        else
            y(i,j) = 0
        endif
    enddo
enddo

return
END SUBROUTINE cx_var

SUBROUTINE cx_var_consec(pc, cp, r, y, whichlin)
type (cx_t_consec), intent (in) :: pc
type (cp_t_consec), dimension(2), intent(in) :: cp
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: y(0:,0:)
integer, intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: i, j
real (kind=dp) :: s,a,c,b, alpha,r0x,d
real (kind=dp), parameter :: pi=3.14159265358979323846_dp
if (size(r,1).ne.size(r,2).or.size(y,1).ne.size(r,1).or. &
  size(y,2).ne.size(r,2)) then
 stop 'cx_var_consec: bad dimensions'
endif

a = r(1,0) !r_HH
b = r(2,0) !r_HC
c = r(2,1) !R_HC

s=(a+b+c)/2._dp
d=dsqrt( (s-a)*(s-b)*(s-c)/s )
alpha= 2._dp*atan( d/(s-a))
if(isnan(alpha)) alpha=pi

do j = 0, size(r,2)-1
    do i = 0, size(r,1)-1
        if (i.ne.j) then
            !! y(i,j) = exp(-r(i,j)/2.6_dp)
            y(i,j) = exp(-r(i,j)/pc%bx)
        else
            y(i,j) = 0
        endif
    enddo
enddo

! if(whichlin.eq.0) then
!     r0x= (-0.269951_dp*alpha+0.435927_dp)*min(b,c) + (0.520510_dp*alpha + 1.737869_dp)
!     y(1,0) = exp(- dabs(r(1,0)-r0x)/pc%bx)*exp(max(r(1,0)-r0x,0._dp)/(pc%cx))
!     y(0,1) = exp(- dabs(r(0,1)-r0x)/pc%bx)*exp(max(r(0,1)-r0x,0._dp)/(pc%cx))
! else if(whichlin.eq.180) then
!     !we have to see, that we take the correct ch distance 
!     !lets just take the hh distance again, because it is unique
!     r0x= (0.8179_dp)*min(b,c) + 4.9980_dp
!     y(1,0) = exp(- dabs(r(1,0)-r0x)/pc%bx)*exp(max(r(1,0)-r0x,0._dp)/(pc%cx))
!     y(0,1) = exp(- dabs(r(0,1)-r0x)/pc%bx)*exp(max(r(0,1)-r0x,0._dp)/(pc%cx))
! else
!     stop 'Illegal whichlin in cx_var_consec'
! end if

if(whichlin.eq.0) then
    r0x= (cp(0)%a * alpha + cp(0)%b) * min(b,c) + (cp(0)%c * alpha + cp(0)%d)
else if(whichlin.eq.180) then
    r0x= (cp(1)%a * alpha + cp(1)%b) * min(b,c) + (cp(1)%c * alpha + cp(1)%d)
else
    stop 'Illegal whichlin in cx_var_consec'
end if

y(1,0) = exp(- dabs(r(1,0)-r0x)/pc%bx)*exp(max(r(1,0)-r0x,0._dp)/(pc%cx))
y(0,1) = exp(- dabs(r(0,1)-r0x)/pc%bx)*exp(max(r(0,1)-r0x,0._dp)/(pc%cx))

return
END SUBROUTINE cx_var_consec

SUBROUTINE cx_write(iun, pc)
! write an element of type cx_t to unit iun
integer, intent (in) :: iun
type (cx_t), intent (in) :: pc

write (iun,'(3(1x,i4),2x,2es20.12,2x,a)') pc, 'type(cx_t)'

return
END SUBROUTINE cx_write

SUBROUTINE cx_write_consec(iun, pc)
! write an element of type cx_t to unit iun
integer, intent (in) :: iun
type (cx_t_consec), intent (in) :: pc

write (iun,'(1(1x,i4),1x,6es20.12,2x,a)') pc, 'type(cx_t_consec)'

return
END SUBROUTINE cx_write_consec

PURE FUNCTION cx_equals(pc0, pc1) RESULT (b)
type (cx_t), intent (in) :: pc0, pc1
logical :: b
!-----------------------------------------------------------------------

b = pc0%dg.eq.pc1%dg.and. &
    pc0%kx.eq.pc1%kx.and. &
    pc0%lx.eq.pc1%lx.and. &
    pc0%ax.eq.pc1%ax.and. &
    pc0%bx.eq.pc1%bx

return
END FUNCTION cx_equals

PURE FUNCTION cx_equals_consec(pc0, pc1) RESULT (b)
type (cx_t_consec), intent (in) :: pc0, pc1
logical :: b
!-----------------------------------------------------------------------

b = pc0%dg.eq.pc1%dg.and. &
    pc0%ax.eq.pc1%ax.and. &
    pc0%kx.eq.pc1%kx.and. &
    pc0%lx.eq.pc1%lx.and. &
    pc0%bx.eq.pc1%bx.and. &
    pc0%cx.eq.pc1%cx.and. &
    pc0%r0x.eq.pc1%r0x

return
END FUNCTION cx_equals_consec

PURE FUNCTION cx_sort(ix) RESULT (iord)
! Stable descending sort, O(n^2) method, for small arrays.
integer, intent (in) :: ix(0:)
integer :: iord(0:size(ix)-1)
!-----------------------------------------------------------------------
integer :: i, k, n, ix0(0:size(ix)-1)

n = size(ix)
ix0 = ix
do i = 0, n-1
    k = maxloc(ix0,dim=1)-1
    iord(i) = k
    ix0(k) = -huge(ix0)
enddo

return
END FUNCTION cx_sort

PURE FUNCTION cx_dim(nb, dg) RESULT (n)
integer, intent (in) :: nb(0:), dg
integer :: n
!-----------------------------------------------------------------------

if (0.le.dg) then
    n = nb(dg)
else
    n = 0
endif

return
END FUNCTION cx_dim

PURE FUNCTION cx_substr(s0, s1) RESULT (b)
! Test if a trimmed version of s1 occurs as a unit in s0.
character (len=*), intent (in) :: s0, s1
logical :: b
!-----------------------------------------------------------------------
integer :: l0, l2
character, parameter :: sep=' '
character (len=len(s1)) :: s2

s2 = adjustl(s1)
if (s0.eq.'*') then
    ! Treat as wildcard
    b = .true.
else
    l0 = len(s0) ; l2 = len_trim(s2)
    if (s0.eq.s2(1:l2)) then
        b = .true.
    else if (l0.le.l2) then
        b = .false.
    else if (s2(1:l2)//sep.eq.s0(1:l2+1).or. &
            sep//s2(1:l2).eq.s0(l0-l2:l0)) then
        b = .true.
    else
        b = index(s0,sep//s2(1:l2)//sep).ne.0
    endif
endif

return
END FUNCTION cx_substr

PURE SUBROUTINE cx_charid(nki, ch)
integer, intent (in) :: nki(0:)
character (len=*), intent (out) :: ch
!-----------------------------------------------------------------------
character, parameter :: &
  digits(0:9)=(/'0','1','2','3','4','5','6','7','8','9'/)
integer :: i, k

if (all(nki.eq.0)) then
    ch = ''
    return
endif
! identify trailing zeros
k = size(nki)
do while (nki(k-1).eq.0)
    k = k-1
enddo
ch = ''
do i = 0, min(len(ch),k)-1
    if (0.le.nki(i).and.nki(i).lt.10) then
        ch(i+1:i+1) = digits(nki(i))
    else
        ch(i+1:i+1) = 'x'
    endif
enddo

return
END SUBROUTINE cx_charid

END MODULE inv_cxx
