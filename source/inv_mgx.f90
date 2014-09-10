MODULE inv_mgx
!..use and access
use inv_core
implicit none

private
public :: mgx_id, mgx_mk1d, mgx_mk2d, mgx_mkrl1d, mgx_mkrl2d, &
  mgx_ngen, mgx_gens, mgx_gens2d, mgx_setd

!..types
type, public :: mgx_t
 integer :: nkk, nk, nr, nkj(0:9), ngrp, ngen, &
   dnpr(0:9), npr(0:9), dnpb(0:9), npb(0:9), dnsc(0:9), nsc(0:9), &
   dnb(0:9), nb(0:9)
 integer, pointer :: gens(:,:), mk1d(:,:)
end type mgx_t

!..procedures
CONTAINS

PURE FUNCTION mgx_id(nkj) RESULT (s)
integer, intent (in) :: nkj(0:)
integer :: s
!-----------------------------------------------------------------------
integer, parameter :: &
  nmax = 12, &
  pt(0:nmax) = (/ 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77 /), &
  ps(0:nmax) = (/ 0, 1, 2, 4, 7, 12, 19, 30, 45, 67, 97, 139, 195 /), &
  pp(0:nmax,0:nmax) = reshape((/ &
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
   0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
   0,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, &
   0,  1,  2,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3, &
   0,  1,  3,  4,  5,  5,  5,  5,  5,  5,  5,  5,  5, &
   0,  1,  3,  5,  6,  7,  7,  7,  7,  7,  7,  7,  7, &
   0,  1,  4,  7,  9, 10, 11, 11, 11, 11, 11, 11, 11, &
   0,  1,  4,  8, 11, 13, 14, 15, 15, 15, 15, 15, 15, &
   0,  1,  5, 10, 15, 18, 20, 21, 22, 22, 22, 22, 22, &
   0,  1,  5, 12, 18, 23, 26, 28, 29, 30, 30, 30, 30, &
   0,  1,  6, 14, 23, 30, 35, 38, 40, 41, 42, 42, 42, &
   0,  1,  6, 16, 27, 37, 44, 49, 52, 54, 55, 56, 56, &
   0,  1,  7, 19, 34, 47, 58, 65, 70, 73, 75, 76, 77/), &
   (/nmax+1,nmax+1/))
integer :: nk, nl, ib, i, j

! Note, pt is the partition function, there "for the record".
! ps are partial sums of the partition function.
! pp(i,j) is the number of partitions of j having length at most i.
nk = sum(nkj)
if (nk.gt.nmax.or.any(nkj.lt.0)) then
    s = -huge(0)
    return
endif
ib = ps(nk)
nl = size(nkj) - count(nkj.eq.0)
i = 1
do while (nk.ne.0)
    ib = ib + pp(nl-1,nk)
    nk = nk - nl
    nl = nl - count(nkj.eq.i)
    i = i + 1
enddo
s = ib

return
END FUNCTION mgx_id

SUBROUTINE mgx_mk1d(nkj, d, x)
! Block revlex order; block sizes nkj(0:)
integer, intent (in) :: nkj(0:)
real (kind=dp), intent (in) :: d(0:,0:)
real (kind=dp), intent (out) :: x(0:)
!-----------------------------------------------------------------------
integer :: i, j, k, l0, l1, j0, j1

if (size(d,1).ne.sum(nkj).or.size(d,2).ne.sum(nkj).or. &
    size(x).ne.sum(nkj)*(sum(nkj)-1)/2) then
    stop 'mgx_mk1d: bad dimensions'
endif

! For reasons of speed we use the code for simple revlex in cases
! where it will give the same result
if (size(nkj).le.1) then
    call mgx_mkrl1d(sum(nkj), d, x)
else if (nkj(1).lt.3.and.all(nkj(2:).lt.2)) then
    call mgx_mkrl1d(sum(nkj), d, x)
else
    k = 0
    j1 = 0
    do l1 = 0, size(nkj)-1
        j0 = 0
        do l0 = 0, l1
            do j = j1, j1+nkj(l1)-1
                do i = j0, min(j-1,j0+nkj(l0)-1)
                    x(k) = d(i,j)
                    k = k + 1
                enddo
            enddo
            j0 = j0 + nkj(l0)
        enddo
        j1 = j1 + nkj(l1)
     enddo
endif

END SUBROUTINE mgx_mk1d

SUBROUTINE mgx_mk2d(nkj, x, d)
! Block revlex order; block sizes nkj(0:)
integer, intent (in) :: nkj(0:)
real (kind=dp), intent (in) :: x(0:)
real (kind=dp) :: d(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j, k, l0, l1, j0, j1

if (size(x).ne.sum(nkj)*(sum(nkj)-1)/2.or. &
    size(d,1).ne.sum(nkj).or.size(d,2).ne.sum(nkj)) then
    stop 'mgx_mk2d: bad dimensions'
endif

! For reasons of speed we use the code for simple revlex in cases
! where it will give the same result
if (size(nkj).le.1) then
    call mgx_mkrl2d(sum(nkj), x, d)
else if (nkj(1).lt.3.and.all(nkj(2:).lt.2)) then
    call mgx_mkrl2d(sum(nkj), x, d)
else
    k = 0
    j1 = 0
    do l1 = 0, size(nkj)-1
        j0 = 0
        do l0 = 0, l1
            do j = j1, j1+nkj(l1)-1
                do i = j0, min(j-1,j0+nkj(l0)-1)
                    d(i,j) = x(k)
                    d(j,i) = d(i,j)
                    k = k + 1
                enddo
            enddo
            j0 = j0 + nkj(l0)
        enddo
        j1 = j1 + nkj(l1)
    enddo
    do i = 0, sum(nkj)-1
        d(i,i) = 0
    enddo
endif

END SUBROUTINE mgx_mk2d

SUBROUTINE mgx_mkrl1d(nk, d, x)
! Simple revlex order
integer, intent (in) :: nk
real (kind=dp), intent (in) :: d(0:,0:)
real (kind=dp), intent (out) :: x(0:)
!-----------------------------------------------------------------------
integer :: i, j, k

if (size(d,1).ne.nk.or.size(d,2).ne.nk.or. &
    size(x).ne.nk*(nk-1)/2) then
    stop 'mgx_mkrl1d: bad dimensions'
endif

k = 0
do j = 0, nk-1
    do i = 0, j-1
        x(k) = d(i,j)
        k = k + 1
    enddo
enddo

END SUBROUTINE mgx_mkrl1d

SUBROUTINE mgx_mkrl2d(nk, x, d)
! Simple revlex order
integer, intent (in) :: nk
real (kind=dp), intent (in) :: x(0:)
real (kind=dp), intent (out) :: d(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j, k

if (size(x).ne.nk*(nk-1)/2.or. &
    size(d,1).ne.nk.or.size(d,2).ne.nk) then
    stop 'mgx_mkrl2d: bad dimensions'
endif

k = 0
do j = 0, nk-1
    do i = 0, j-1
        d(i,j) = x(k)
        d(j,i) = d(i,j)
        k = k + 1
    enddo
enddo
do i = 0, nk-1
    d(i,i) = 0
enddo

END SUBROUTINE mgx_mkrl2d

PURE FUNCTION mgx_ngen(nkj) RESULT (ng)
integer, intent (in) :: nkj(0:)
integer :: ng
!-----------------------------------------------------------------------
integer :: i, n0

n0 = 0
do i = 0, size(nkj)-1
    if (3.le.nkj(i)) then
        n0 = n0+2
    else if (nkj(i).eq.2) then
        n0 = n0+1
    endif
enddo
ng = n0

return
END FUNCTION mgx_ngen

SUBROUTINE mgx_gens(nkj, ind, iord)
integer, intent (in) :: nkj(0:), ind
integer, intent (out) :: iord(0:)
!-----------------------------------------------------------------------
integer :: nk, i
integer, allocatable :: iord0(:)
real (kind=dp), allocatable :: x(:), d(:,:)

nk = sum(nkj)
allocate (iord0(0:nk-1))
allocate (x(0:nk*(nk-1)/2-1), d(0:nk-1,0:nk-1))
if (size(iord).ne.nk*(nk-1)/2) then
    stop 'mgx_gens: bad size iord'
endif
call mgx_gens2d(nkj, ind, iord0)
do i = 0, nk*(nk-1)/2-1
    x(i) = i
enddo
call mgx_mk2d(nkj, x, d)
call mgx_mk1d(nkj, d(iord0,iord0), x)
do i = 0, nk*(nk-1)/2-1
    iord(i) = x(i)
enddo

return
END SUBROUTINE mgx_gens

SUBROUTINE mgx_gens2d(nkj, ind, iord)
integer, intent (in) :: nkj(0:), ind
integer, intent (out) :: iord(0:)
!-----------------------------------------------------------------------
integer :: i, k0, n0, k

if (size(iord).ne.sum(nkj)) then
    stop 'mgx_gens2d: bad size iord'
endif

n0 = 0; k0 = 0
do i = 0, size(nkj)-1
    if (nkj(i).lt.2.or. &
        (ind.ne.n0.and.(nkj(i).eq.2.or.ind.ne.n0+1))) then
        iord(k0:k0+nkj(i)-1) = (/ (k0+k, k=0,nkj(i)-1) /)
    else if (ind.eq.n0) then
        iord(k0:k0+nkj(i)-1) = (/ k0+1, k0, (k0+k, k=2,nkj(i)-1) /)
    else if (ind.eq.n0+1) then
        iord(k0:k0+nkj(i)-1) = (/ (k0+modulo(k,nkj(i)), k=1,nkj(i)) /)
    endif
    if (3.le.nkj(i)) then
        n0 = n0 + 2
    else if (nkj(i).eq.2) then
        n0 = n0 + 1
    endif
    k0 = k0 + nkj(i)
enddo

return
END SUBROUTINE mgx_gens2d

SUBROUTINE mgx_setd(r, d, d2, d3, d4, d5, d6, d7, d8, d9)
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: d(0:,0:)
real (kind=dp), optional, intent (out) :: d2(0:,0:), d3(0:,0:), &
  d4(0:,0:), d5(0:,0:), d6(0:,0:), d7(0:,0:), d8(0:,0:), d9(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j

if (size(d,1).ne.size(r,1).or.size(d,2).ne.size(r,2).or. &
    size(r,1).ne.size(r,2)) then
    stop 'mgx_setd: bad dimensions'
endif

do j = 0, size(r,2)-1
    do i = 0, j-1
        d(i,j) = (r(i,j) + r(j,i)) / 2._dp
        d(j,i) = d(i,j)
    enddo
        d(j,j) = 0._dp
enddo
if (present(d2)) then
    d2 = d * d
    if (present(d3)) then
        d3 = d2 * d
        if (present(d4)) then
            d4 = d3 * d
            if (present(d5)) then
                d5 = d4 * d
                if (present(d6)) then
                    d6 = d5 * d
                    if (present(d7)) then
                        d7 = d6 * d
                        if (present(d8)) then
                            d8 = d7 * d
                            if (present(d9)) then
                                d9 = d8 * d
                            endif
                        endif
                    endif
                endif
            endif
        endif
    endif
endif

return
END SUBROUTINE mgx_setd

END MODULE inv_mgx
