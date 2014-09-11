PURE FUNCTION px_x2_nbase () RESULT (nb)
integer :: nb
!-----------------------------------------------------------------------
nb = cx_nbase(pes_x2_nki,pes_x2_sysnew,px_pcv)
return
END FUNCTION px_x2_nbase
