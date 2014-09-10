MODULE inv_cxz
!..use and access
use inv_dp
use inv_core
use inv_mg
use inv_cxx
use inv_cxy
implicit none

private
public :: cx_nbx, cx_bx, cx_fx, cx_nbase, cx_base, cx_getcf
public :: cx_bx_consec, cx_nbase_consec, cx_base_consec, cx_getcf_consec

!..procedures
CONTAINS

PURE FUNCTION cx_nbx(nkj, dg) RESULT (nb)
integer, intent (in) :: nkj(0:), dg
integer :: nb
!-----------------------------------------------------------------------
! Note: All error conditions are indicated by returning nb = -huge(0)
integer :: id

if (dg.lt.0) then
    nb = 0
    return
endif
select case (mgx_id(nkj))
    case (0)
        nb = 0
    case (1)
        nb = 1 !! must look into this
    case (mg2_id)
        nb = mg2_nb(dg)
    case (mg11_id)
        nb = mg11_nb(dg)
    case (mg21_id)
        nb = mg21_nb(dg)
    case default
        nb = -huge(0)
end select

return
END FUNCTION cx_nbx

SUBROUTINE cx_bx(nki, nkj, pc, r, w)
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
integer :: ikv(0:size(nkj)-1), len

if (size(nki).lt.size(nkj)) then
    stop 'cx_bx: bad dimension nki, nkj'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
    stop 'cx_bx: bad dimension nki, r'
!! else if (size(w).ne....) then
!!  stop 'cx_bx: bad dimension w'
endif

if (any(nki(0:size(nkj)-1).lt.nkj)) then
    w = 0
    return
else if (2.le.sum(nkj).and.pc%dg.lt.0) then
    w = 0
    return
endif

ikv = cx_sort(nkj)
len = count(nkj.ne.0)
w = 0
select case (mgx_id(nkj))
    case (0)
    case (1)
        call cx_b1 (nki, ikv(0:len-1), w)
    case (mg2_id)
        call cx_b2 (nki, ikv(0:len-1), pc, r, w)
    case (mg11_id)
        call cx_b11 (nki, ikv(0:len-1), pc, r, w)
    case (mg21_id)
        call cx_b21 (nki, ikv(0:len-1), pc, r, w)
    case default
        stop 'cx_bx: bad nkj'
end select

return
END SUBROUTINE cx_bx

SUBROUTINE cx_bx_consec(nki, nkj, pc, cp, r, w, whichlin)
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t_consec), intent (in) :: pc
type (cp_t_consec), dimension(2), intent(in) :: cp
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:)
integer, intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: ikv(0:size(nkj)-1), len

if (size(nki).lt.size(nkj)) then
    stop 'cx_bx_consec: bad dimension nki, nkj'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
    stop 'cx_bx_consec: bad dimension nki, r'
!! else if (size(w).ne....) then
!!  stop 'cx_bx: bad dimension w'
endif

if (any(nki(0:size(nkj)-1).lt.nkj)) then
    w = 0
    return
else if (2.le.sum(nkj).and.pc%dg.lt.0) then
    w = 0
    return
endif

ikv = cx_sort(nkj)
len = count(nkj.ne.0)
w = 0
select case (mgx_id(nkj))
    !case (0)
    case (mg21_id)
        call cx_b21_consec (nki, ikv(0:len-1), pc, cp, r, w,whichlin)
    case default
        stop 'cx_bx_consec: bad nkj'
end select

return
END SUBROUTINE cx_bx_consec

FUNCTION cx_fx(nki, nkj, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:), nkj(0:)
type (cx_t), intent (in) :: pc
real (kind=dp), intent (in) :: r(0:,0:), cf(0:)
real (kind=dp) :: f
!-----------------------------------------------------------------------
real (kind=dp) :: w(0:size(cf)-1)

call cx_bx(nki, nkj, pc, r, w)
f = dot_product(cf, w)

return
END FUNCTION cx_fx

PURE FUNCTION cx_nbase(nki, sys, pcv) RESULT (nb)
integer, intent (in) :: nki(0:)
character (len=*), intent (in) :: sys
type (cx_t), intent (in) :: pcv(2:)
integer :: nb
!-----------------------------------------------------------------------
integer :: nkj(0:size(nki)-1), order, ib, i, k, l
logical :: done
character (len=size(nki)) :: charid

! initialize counter ib
ib = 0
! One-body terms
do i = 0, size(nki)-1
    if (1.le.nki(i)) then
        nkj = 0 ; nkj(i) = 1
        call cx_charid(nkj, charid)
        if (cx_substr(sys,charid)) then
            ib = ib + 1
        endif
    endif
enddo

! Further terms grevlex
do order = 2, min(sum(nki),ubound(pcv,dim=1))
    if (0.le.pcv(order)%dg) then
        k = order
        do i = 0, size(nki)-1
            nkj(i) = min(k, nki(i))
            k = k - nkj(i)
        enddo
        done = .false.
        do while (.not.done)
            call cx_charid(nkj, charid)
            if (cx_substr(sys,charid)) then
                ib = ib + cx_nbx(nkj,pcv(order)%dg)
            endif
            ! set next nkj
            l = size(nki)
            do i = size(nki)-1, 0, -1
                if (nkj(i).lt.nki(i).and.sum(nkj(i:)).lt.order) then
                    l = i
                endif
            enddo
            if (l.lt.size(nki)) then
                nkj(l) = nkj(l)+1
                k = order-sum(nkj(l:))
                do i = 0, l-1
                    nkj(i) = min(k, nki(i))
                    k = k - nkj(i)
                enddo
                done = .false.
            else
                done = .true.
            endif
        enddo
    endif
enddo
nb = ib

return
END FUNCTION cx_nbase

PURE FUNCTION cx_nbase_consec(nki, sys, pcv) RESULT (nb)
integer, intent (in) :: nki(0:)
character (len=*), intent (in) :: sys
type (cx_t_consec), intent (in) :: pcv(2:)
integer :: nb
!-----------------------------------------------------------------------
integer :: nkj(0:size(nki)-1), order, ib, i, k, l
logical :: done
character (len=size(nki)) :: charid

! initialize counter ib
ib = 0
! One-body terms
do i = 0, size(nki)-1
    if (1.le.nki(i)) then
        nkj = 0 ; nkj(i) = 1
        call cx_charid(nkj, charid)
        if (cx_substr(sys,charid)) then
            ib = ib + 1
        endif
    endif
enddo

! Further terms grevlex
do order = 2, min(sum(nki),ubound(pcv,dim=1))
    if (0.le.pcv(order)%dg) then
        k = order
        do i = 0, size(nki)-1
            nkj(i) = min(k, nki(i))
            k = k - nkj(i)
        enddo
        done = .false.
        do while (.not.done)
            call cx_charid(nkj, charid)
            if (cx_substr(sys,charid)) then
                ib = ib + cx_nbx(nkj,pcv(order)%dg)
            endif
            ! set next nkj
            l = size(nki)
            do i = size(nki)-1, 0, -1
                if (nkj(i).lt.nki(i).and.sum(nkj(i:)).lt.order) then
                    l = i
                endif
            enddo
            if (l.lt.size(nki)) then
                nkj(l) = nkj(l)+1
                k = order-sum(nkj(l:))
                do i = 0, l-1
                    nkj(i) = min(k, nki(i))
                    k = k - nkj(i)
                enddo
                done = .false.
            else
                done = .true.
            endif
        enddo
    endif
enddo
nb = ib

return
END FUNCTION cx_nbase_consec

SUBROUTINE cx_base(nki, sys, pcv, r, w)
integer, intent (in) :: nki(0:)
character (len=*), intent (in) :: sys
type (cx_t), intent (in) :: pcv(2:)
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
integer :: nkj(0:size(nki)-1), order, ib, nb, i, k, l
logical :: done
character (len=size(nki)) :: charid

if (size(w).ne.cx_nbase(nki,sys,pcv)) then
    stop 'cx_base: bad dimension w'
endif

! initialize counter ib
ib = 0
! One-body terms
do i = 0, size(nki)-1
    if (1.le.nki(i)) then
        nkj = 0 ; nkj(i) = 1
        call cx_charid (nkj, charid)
        if (cx_substr(sys,charid)) then
            w(ib) = nki(i)
            ib = ib + 1
        endif
    endif
enddo

! Further terms grevlex
do order = 2, min(sum(nki),ubound(pcv,dim=1))
    if (0.le.pcv(order)%dg) then
        k = order
        do i = 0, size(nki)-1
            nkj(i) = min(k, nki(i))
            k = k - nkj(i)
        enddo
        done = .false.
        do while (.not.done)
            call cx_charid(nkj, charid)
            if (cx_substr(sys,charid)) then
                nb = cx_nbx(nkj,pcv(order)%dg)
                call cx_bx(nki, nkj, pcv(order), r, w(ib:ib+nb-1))
                ib = ib + nb
            endif
            ! set next nkj
            l = size(nki)
            do i = size(nki)-1, 0, -1
                if (nkj(i).lt.nki(i).and.sum(nkj(i:)).lt.order) then
                    l = i
                endif
            enddo
            if (l.lt.size(nki)) then
                nkj(l) = nkj(l) + 1
                k = order - sum(nkj(l:))
                do i = 0, l-1
                    nkj(i) = min(k, nki(i))
                    k = k - nkj(i)
                enddo
                done = .false.
            else
                done = .true.
            endif
        enddo
    endif
enddo

if (ib.ne.size(w)) then
    stop 'cx_base: size mismatch'
endif

return
END SUBROUTINE cx_base

SUBROUTINE cx_base_consec(nki, sys, pcv, cp, r, w, whichlin)
integer, intent (in) :: nki(0:)
character (len=*), intent (in) :: sys
type (cx_t_consec), intent (in) :: pcv(2:)
type (cp_t_consec), dimension(2), intent(in) :: cp
real (kind=dp), intent (in) :: r(0:,0:)
real (kind=Dp), intent (out) :: w(0:)
integer,intent(in) :: whichlin
!-----------------------------------------------------------------------
integer :: nkj(0:size(nki)-1), order, ib, nb, i, k, l
logical :: done
character (len=size(nki)) :: charid

if (size(w).ne.cx_nbase_consec(nki,sys,pcv)) then
    stop 'cx_base: bad dimension w'
endif

! initialize counter ib
ib = 0
! One-body terms
do i = 0, size(nki)-1
    if (1.le.nki(i)) then
        nkj = 0 ; nkj(i) = 1
        call cx_charid(nkj, charid)
        if (cx_substr(sys,charid)) then
            w(ib) = nki(i)
            ib = ib + 1
        endif
    endif
enddo
! Further terms grevlex
do order = 2, min(sum(nki),ubound(pcv,dim=1))
    if (0.le.pcv(order)%dg) then
        k = order
        do i = 0, size(nki)-1
            nkj(i) = min(k,nki(i))
            k = k - nkj(i)
        enddo
        done = .false.
        do while (.not.done)
            call cx_charid(nkj, charid)
            if (cx_substr(sys,charid)) then
                nb = cx_nbx(nkj,pcv(order)%dg)
                call cx_bx_consec(nki, nkj, pcv(order), cp, r, w(ib:ib+nb-1), whichlin)
                ib = ib + nb
            endif
            ! set next nkj
            l = size(nki)
            do i = size(nki)-1, 0, -1
                if (nkj(i).lt.nki(i).and.sum(nkj(i:)).lt.order) then
                    l = i
                endif
            enddo
            if (l.lt.size(nki)) then
                nkj(l) = nkj(l) + 1
                k = order - sum(nkj(l:))
                do i = 0, l-1
                    nkj(i) = min(k, nki(i))
                    k = k - nkj(i)
                enddo
                done = .false.
            else
                done = .true.
            endif
        enddo
    endif
enddo

if (ib.ne.size(w)) then
    stop 'cx_base_consec: size mismatch'
endif

return
END SUBROUTINE cx_base_consec

SUBROUTINE cx_getcf(proc, nki, sys, pcv, cf)
! Extract the separate coefficient blocks from array cf.
interface
    subroutine proc(nki, nkj, pc, cf)
    use inv_dp
    use inv_cxx
    integer, intent (in) :: nki(0:), nkj(0:)
    type (cx_t), intent (in) :: pc
    real (kind=dp), intent (in) :: cf(0:)
    end subroutine proc
end interface
integer, intent (in) :: nki(0:)
character (len=*), intent (in) :: sys
type (cx_t), intent (in) :: pcv(2:)
real (kind=dp), intent (in) :: cf(0:)
!-----------------------------------------------------------------------
integer :: l0, nkj(0:size(nki)-1), order, ib, nb, i, k, l
logical :: done
character (len=size(nki)) :: charid

if (size(cf).ne.cx_nbase(nki,sys,pcv)) then
    stop 'cx_getcf: bad dimension cf'
endif

l0 = len_trim(sys)
! initialize counter ib
ib = 0
! One-body terms
do i = 0, size(nki)-1
    if (1.le.nki(i)) then
        nkj = 0 ; nkj(i) = 1
        call cx_charid(nkj, charid)
        if (cx_substr(sys(1:l0),charid)) then
            nb = 1
            call proc(nki, nkj, cx_null, cf(ib:ib+nb-1))
            ib = ib + nb
        endif
    endif
enddo
! Further terms grevlex
do order = 2, min(sum(nki),ubound(pcv,dim=1))
    k = order
    do i = 0, size(nki)-1
        nkj(i) = min(k,nki(i))
        k = k - nkj(i)
    enddo
    done = .false.
    do while (.not.done)
        call cx_charid(nkj, charid)
        if (cx_substr(sys(1:l0),charid)) then
            nb = cx_nbx(nkj,pcv(order)%dg)
            call proc(nki, nkj, pcv(order), cf(ib:ib+nb-1))
            ib = ib + nb
        endif
        ! set next nkj
        l = size(nki)
        do i = size(nki)-1, 0, -1
            if (nkj(i).lt.nki(i).and.sum(nkj(i:)).lt.order) then
                l = i
            endif
        enddo
        if (l.lt.size(nki)) then
            nkj(l) = nkj(l) + 1
            k = order - sum(nkj(l:))
            do i = 0, l-1
                nkj(i) = min(k,nki(i))
                k = k - nkj(i)
            enddo
            done = .false.
        else
            done = .true.
        endif
    enddo
enddo

if (ib.ne.size(cf)) then
    stop 'cx_getcf: size mismatch'
endif

return
END SUBROUTINE cx_getcf

SUBROUTINE cx_getcf_consec(proc, nki, sys, pcv, cf, whichlin)
! Extract the separate coefficient blocks from array cf.
interface
    subroutine proc (nki, nkj, pc, cf,whichlin)
    use inv_dp
    use inv_cxx
    integer, intent (in) :: nki(0:), nkj(0:)
    type (cx_t_consec), intent (in) :: pc
    real (kind=dp), intent (in) :: cf(0:)
    integer, intent (in) :: whichlin
    end subroutine proc
end interface
integer, intent (in) :: nki(0:)
character (len=*), intent (in) :: sys
type (cx_t_consec), intent (in) :: pcv(2:)
real (kind=dp), intent (in) :: cf(0:)
integer, intent (in) :: whichlin
!-----------------------------------------------------------------------
integer :: l0, nkj(0:size(nki)-1), order, ib, nb, i, k, l
logical :: done
character (len=size(nki)) :: charid

if (size(cf).ne.cx_nbase_consec(nki,sys,pcv)) then
    stop 'cx_getcf_consec: bad dimension cf'
endif

l0 = len_trim(sys)
! initialize counter ib
ib = 0
! One-body terms
do i = 0, size(nki)-1
    if (1.le.nki(i)) then
        nkj = 0 ; nkj(i) = 1
        call cx_charid(nkj, charid)
        if (cx_substr(sys(1:l0),charid)) then
            nb = 1
            call proc(nki, nkj, cx_null_consec, cf(ib:ib+nb-1), whichlin)
            ib = ib + nb
        endif
    endif
enddo
! Further terms grevlex
do order = 2, min(sum(nki),ubound(pcv,dim=1))
    k = order
    do i = 0, size(nki)-1
        nkj(i) = min(k,nki(i))
        k = k - nkj(i)
    enddo
    done = .false.
    do while (.not.done)
        call cx_charid(nkj, charid)
        if (cx_substr(sys(1:l0),charid)) then
            nb = cx_nbx(nkj,pcv(order)%dg)
            call proc(nki, nkj, pcv(order), cf(ib:ib+nb-1), whichlin)
            ib = ib + nb
        endif
        ! set next nkj
        l = size(nki)
        do i = size(nki)-1, 0, -1
            if (nkj(i).lt.nki(i).and.sum(nkj(i:)).lt.order) then
                l = i
            endif
        enddo
        if (l.lt.size(nki)) then
            nkj(l) = nkj(l)+1
            k = order-sum(nkj(l:))
            do i = 0, l-1
                nkj(i) = min(k, nki(i))
                k = k - nkj(i)
            enddo
            done = .false.
        else
            done = .true.
        endif
    enddo
enddo

if (ib.ne.size(cf)) then
    stop 'cx_getcf_consec: size mismatch'
endif

return
END SUBROUTINE cx_getcf_consec

END MODULE inv_cxz
