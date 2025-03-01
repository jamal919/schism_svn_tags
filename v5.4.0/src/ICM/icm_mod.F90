!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

module icm_mod
!-------------------------------------------------------------------------------
!parameter definition for ICM
!-------------------------------------------------------------------------------
  use schism_glbl,only: rkind,nea,nvrt
  implicit none

  !integer, parameter ::iT=2
  !real(kind=rkind), parameter :: CV1=1.0E8
  !real(kind=rkind), parameter :: CV2=1.0E8
  real(kind=rkind), parameter :: COV=1.0E-10
  !molar weight for C,Ca,CaCo3,N
  real(kind=rkind), parameter :: mC=12.011,mCACO3=100.086,mN=14.007


  !time step in ICM
  real(kind=rkind), save :: dtw,dtw2

  !time stamp for WQinput
  real(kind=rkind),save:: time_icm(5),time_ph
  
  !global switch 
  integer,save :: iSun,iNPS,iPS,iLight,jLight,iSed,iRea,iZoo,iPh
  integer,save :: iAtm,iBen,iRad,iCheck,iout_icm
  integer,save :: iWS,iTurb,iWRea,iTSS  

  !water quality state variables
  real(kind=rkind),save,allocatable,dimension(:,:,:) :: wqc
  real(kind=rkind),save,allocatable,dimension(:) :: dep,Sal,Temp,TSED
  real(kind=rkind),save,allocatable,dimension(:,:) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOCA,RPON,LPON,DON,NH4,NO3
  real(kind=rkind),save,allocatable,dimension(:,:) :: RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOC

  !PH model
  integer, save :: inu_ph,irec_ph
  integer,save,allocatable :: iphgb(:)
  real(kind=rkind),save,allocatable :: ph_nudge(:),ph_nudge_nd(:) 
  real(kind=rkind),save,allocatable,dimension(:,:) :: TIC,ALK,CA,CACO3,PH_el,PH_nd,TIC_el,ALK_el                         
  real(kind=rkind),save,allocatable,dimension(:) :: PH,CAsat,CO2
  real(kind=rkind),save :: WSCACO3,rKCACO3,rKCA,rKa

  !phyto. growth rate
  real(kind=rkind),save :: TU,TD,rIa,Daylen,rIavg
  real(kind=rkind),save,allocatable,dimension(:,:) :: GP,PrefN

  !TSED
  real(kind=rkind),save,allocatable,dimension(:) :: PC2TSS 
  real(kind=rkind),save :: WSSED
  
  !DO
  real(kind=rkind),save,allocatable,dimension(:) :: WMS 

  !---------general parameters from icm.in--------------------------------
  !zooplankton paramters
  real(kind=rkind),save :: Eff,RF,Pf
  real(kind=rkind),save,dimension(8,2) :: GZM,rKhGE,PPC
  real(kind=rkind),save,dimension(2) :: BMZR,DRZ,TGZ,rKTGZ1,rKTGZ2,TBZ,rKTBZ,RZ

  !phytoplankton parameters 
  real(kind=rkind),save :: rKhS,ST,rKeC1,rKeC2,mKhN,mKhP,Dopt 
  real(kind=rkind),save,dimension(3) :: GPM,BMPR,PRR,TGP,rKTGP1,rKTGP2,TBP,rKTBP,CChl,rKhN,rKhP,rIm,alpha_PB

  !carbon parameters 
  real(kind=rkind),save :: FCRPZ,FCLPZ,FCDPZ,FCRP,FCLP,FCDP
  real(kind=rkind),save :: rKRC,rKLC,rKDC,rKRCalg,rKLCalg,rKDCalg,TRHDR,TRMNL,rKTHDR,rKTMNL
  real(kind=rkind),save :: rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
  real(kind=rkind),save,dimension(3) :: FCD
  real(kind=rkind),save,dimension(2) :: FCDZ,rKHRZ

  !nitrogen parameters 
  real(kind=rkind),save :: FNRPZ,FNLPZ,FNDPZ,FNIPZ,FNRP,FNLP,FNDP,FNIP,ANDC
  real(kind=rkind),save :: rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg,rNitM,TNit,rKNit1,rKNit2,rKhNitDO,rKhNitN
  real(kind=rkind),save,dimension(3) :: FNR,FNL,FND,FNI,ANC
  real(kind=rkind),save,dimension(2) :: FNRZ,FNLZ,FNDZ,FNIZ,ANCZ

  !phosphorus parameters 
  real(kind=rkind),save :: FPRPZ,FPLPZ,FPDPZ,FPIPZ,FPRP,FPLP,FPDP,FPIP
  real(kind=rkind),save :: rKPO4p,rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg 
  real(kind=rkind),save,dimension(3) :: FPR,FPL,FPD,FPI,APC
  real(kind=rkind),save,dimension(2) :: FPRZ,FPLZ,FPDZ,FPIZ,APCZ

  !silica parameters 
  real(kind=rkind),save :: FSPPZ,FSIPZ,FSPP,FSIP,rKSAp,rKSU,TRSUA,rKTSUA
  real(kind=rkind),save :: FSPd,FSId,ASCd
  real(kind=rkind),save,dimension(2) :: FSPZ,FSIZ,ASCZ

  !COD&DO parameters 
  real(kind=rkind),save :: rKHCOD,rKCD,TRCOD,rKTCOD  
  real(kind=rkind),save :: AOC,AON,rKro,rKTr         
  !--------------------------------------------------------------------------------------

  !settling
  real(kind=rkind),save,allocatable,dimension(:) :: WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea

  !benthic flux from sediment flux model
  real(kind=rkind),save:: BnDOC,BnNH4,BnNO3,BnPO4t,BnSAt,BnCOD,BnDO

  !additional benthic flux 
  real(kind=rkind),save:: TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
  real(kind=rkind),save,allocatable,dimension(:) :: BRPOC,BLPOC,BDOC,BRPON,BLPON,BDON,BNH4,BNO3,BRPOP,BLPOP,BDOP,BPO4t,BSU,BSAt,BCOD,BDO

  !surface flux : atmospheric loading
  real(kind=rkind),save :: SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3,SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO

  !loading
  real(kind=rkind),save,allocatable,dimension(:) :: WWPRPOC,WWPLPOC,WWPDOC,WWPRPON,WWPLPON,WWPDON,WWPNH4,WWPNO3,&
                                                   & WWPRPOP,WWPLPOP,WWPDOP,WWPPO4t,WWPSU,WWPSAt,WWPCOD,WWPDO,WWPSalt 
  real(kind=rkind),save:: WPRPOC,WPLPOC,WPDOC,WPRPON,WPLPON,WPDON,WPNH4,WPNO3,WPRPOP,WPLPOP,WPDOP,WPPO4t,WPSU,WPSAt,WPCOD,WPDO 
  real(kind=rkind),save :: WZB1,WZB2,WPB1,WPB2,WPB3,WRPOC,WLPOC,WDOC,WRPON,WLPON,WDON,WNH4,WNO3,WRPOP,WLPOP,WDOP,WPO4t,WSU,WSAt,WCOD,WDO      

  !for station output for intermediate parameters and ICM variables
  !ista(ie) refers to local station index (lsi)
  !nsta(lsi) refers to number of depth
  !depsta(k,lsi) is depth,where k is depth index
  !stanum is the station index from cstation.in
  integer, save :: nspool_icm
  integer,save,allocatable :: ista(:),nsta(:),stanum(:,:)
  real(rkind),save,allocatable :: depsta(:,:)

end module icm_mod
