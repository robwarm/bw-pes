SUBROUTINE px_null_base (xn0, w)
real (kind=dp), intent (in) :: xn0(0:,0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
if (size(w).ne.px_null_nbase()) then
 stop 'px_null_base: bad dimension'
endif
w = 0
return
END SUBROUTINE px_null_base
