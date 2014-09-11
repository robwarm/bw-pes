SUBROUTINE px_x2y1consec_base (xn0, w)
real (kind=dp), intent (in) :: xn0(0:,0:)
real (kind=dp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
real (kind=dp) :: r(0:size(xn0,2)-1,0:size(xn0,2)-1)

real(kind=dp) :: nb_consec0,nb_consec180, nb_tot

if (size(w).ne.(px_x2y1_nbase()+px_x2y1consec_nbase(px_pcv_consec0)+px_x2y1consec_nbase(px_pcv_consec180))) then
 stop 'px_x2y1consec_base: bad dimension'
endif

nb_tot=size(w)
nb_consec0= cx_nbase_consec(pes_x2y1consec_nki,pes_x2y1consec_sysnew,px_pcv_consec0)
nb_consec180= cx_nbase_consec(pes_x2y1consec_nki,pes_x2y1consec_sysnew,px_pcv_consec180)

call pes_dists (xn0, r)
call cx_base (pes_x2y1_nki, pes_x2y1_sysnew, px_pcv, r, w(0:(nb_tot-1-nb_consec0-nb_consec180)))

call cx_base_consec (pes_x2y1consec_nki, pes_x2y1consec_sysnew, px_pcv_consec0, pes_x2y1consec_cp, r, w(nb_tot-nb_consec0-nb_consec180:nb_tot-1-nb_consec180),0)
call cx_base_consec (pes_x2y1consec_nki, pes_x2y1consec_sysnew, px_pcv_consec180, pes_x2y1consec_cp, r, w(nb_tot-nb_consec180:nb_tot-1),180)
return
END SUBROUTINE px_x2y1consec_base
