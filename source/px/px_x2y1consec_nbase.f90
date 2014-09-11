PURE FUNCTION px_x2y1consec_nbase (px_pcv_c) RESULT (nb)
type (cx_t_consec), intent(in) :: px_pcv_c(2:)
integer :: nb
!-----------------------------------------------------------------------
nb =  cx_nbase_consec(pes_x2y1consec_nki,pes_x2y1consec_sysnew,px_pcv_c)
return
END FUNCTION px_x2y1consec_nbase
