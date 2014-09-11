MODULE px
!..use and access
use inv
use pes
implicit none
private
public :: &
  px_lsq, px_lsqs, px_errf, &
  px_null_nbase, px_null_base
public &
  px_x2_nbase, px_x2_base, &
  px_x1y1_nbase, px_x1y1_base, &
  px_x2y1_nbase, px_x2y1_base

public &
  px_x2y1consec_nbase, px_x2y1consec_base

!..data
save
character (len=*), parameter, public :: &
  pes_x2_sysold_dflt='x1', &
  pes_x2_sysnew_dflt='2', &
  pes_x1y1_sysold_dflt='x1 y1', &
  pes_x1y1_sysnew_dflt='11', &
  pes_x2y1_sysold_dflt='x1 y1', &
  pes_x2y1_sysnew_dflt='2 11 21'

character (len=*), parameter, public :: &
  pes_x2y1consec_sysold_dflt='', &
  pes_x2y1consec_sysnew_dflt='21'

character (len=240), public :: &
  pes_x2_sysold=pes_x2_sysold_dflt, &
  pes_x2_sysnew=pes_x2_sysnew_dflt, &
  pes_x1y1_sysold=pes_x1y1_sysold_dflt, &
  pes_x1y1_sysnew=pes_x1y1_sysnew_dflt, &
  pes_x2y1_sysold=pes_x2y1_sysold_dflt, &
  pes_x2y1_sysnew=pes_x2y1_sysnew_dflt

character (len=240), public :: &
  pes_x2y1consec_sysold=pes_x2y1consec_sysold_dflt, &
  pes_x2y1consec_sysnew=pes_x2y1consec_sysnew_dflt

integer, public :: px_ng=3, px_ng2=3
real (kind=dp), public :: &
  px_rcond=1.0e-11_dp, &
  px_wtgf=0.2_dp*pes_bohr, &
  px_wtg2f=0.05_dp*pes_bohr**2
type (cx_t), public :: px_pcv(2:3)= &
  (/ cx_null, cx_null /)

type (cx_t_consec), public :: px_pcv_consec0(2:3)= &
  (/ cx_null_consec, cx_null_consec /)

type (cx_t_consec), public :: px_pcv_consec180(2:3)= &
  (/ cx_null_consec, cx_null_consec /)

!..procedures
CONTAINS

include 'px_lsq.f90'
include 'px_lsqs.f90'
include 'px_errf.f90'
include 'px_null_nbase.f90'
include 'px_null_base.f90'
include 'px_x2_nbase.f90'
include 'px_x2_base.f90'
include 'px_x1y1_nbase.f90'
include 'px_x1y1_base.f90'
include 'px_x2y1_nbase.f90'
include 'px_x2y1_base.f90'

include 'px_x2y1consec_nbase.f90'
include 'px_x2y1consec_base.f90'

END MODULE px
