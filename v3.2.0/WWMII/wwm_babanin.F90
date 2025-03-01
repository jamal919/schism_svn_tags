MODULE SdsBabanin

CONTAINS

  subroutine calc_Sds(ip,nfreq,EDENS,f,Kds,ANAR_IN,LPOW,MPOW,a1,a2,ACLOC,IMATRA,IMATDA,SSDS)
    USE DATAPOOL
    IMPLICIT NONE

    integer, intent(in)            :: ip

    real             , INTENT(IN)  ::  EDENS(:)   ! E(f) m2/Hz
    real             , INTENT(IN)  ::  f(:)       ! f in Hz
    real             , INTENT(IN)  ::  ANAR_IN(:) ! directional narrowness as defined in Babanin publications (input value)
    real             , INTENT(IN)  ::  a1         ! coefficient on T1
    real             , INTENT(IN)  ::  a2         ! coefficient on T2
    real             , INTENT(IN)  ::  LPOW       ! power on threshold exceedence in T1 (precise definition depends on whether using "U" vs "D" method)
    real             , INTENT(IN)  ::  MPOW       ! power on threshold exceedence in T2 (precise definition depends on whether using "U" vs "D" method)
    integer          , INTENT(IN)  ::  nfreq      ! # freqs
    real             , INTENT(OUT) ::  SSDS(MSC,MDC) !
    real             , INTENT(IN)  :: ACLOC(MSC,MDC)
    real             , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
    real             , INTENT(OUT) ::  Kds(:)     ! Kds(f)=Sds(f)/E(f)


    real             , allocatable ::  Sds(:)     ! Sds(f), the source term
    real             , allocatable ::  ndedens(:) ! NDEDENS(f)=DEDENS(f)/EDENST(f)
    real             , allocatable ::  dedens(:)  ! DEDENS(f)=EDENS(f)-EDENST(f)
    real             , allocatable ::  edenst(:)  ! EDENST(f)=(2*g.^2)*((2*pi).^-4)*(f.^-5).*(A.^-1)*Bnt
    real             , allocatable ::  T1(:)      ! inherent dissipation (divided by E(f))
    real             , allocatable ::  T2(:)      ! induced dissipation (divided by E(f))
    real             , allocatable ::  ST1(:)     ! inherent dissipation
    real             , allocatable ::  ST2(:)     ! induced dissipation
    real             , allocatable ::  ANAR(:)    ! directional narrowness as defined in Babanin publications (value used)
    real             , allocatable ::  xff(:),ADF(:) ! temporary arrays for integration

    REAL ::  ASUM   ! temporary variable for integration
    REAL ::  BNT    ! is an empirical coefficient related to the spectral density in the saturation range (Banner et al 2007)
    real ::  ST1_INT,ST2_INT,Sds_int ! integrated values (for test output only)
    INTEGER :: II,IS,ID ! counters
!   real ::  ELIM ! needed for 42D (can comment but don't delete)
    
    ALLOCATE(sds(NFREQ))
    ALLOCATE(xff(NFREQ))
    ALLOCATE(ndedens(NFREQ))
    ALLOCATE(dedens(NFREQ))
    ALLOCATE(edenst(NFREQ))
    ALLOCATE(T1(NFREQ))
    ALLOCATE(T2(NFREQ))
    ALLOCATE(ANAR(NFREQ))
    ALLOCATE(ST1(NFREQ))
    ALLOCATE(ST2(NFREQ))
    ALLOCATE(ADF(NFREQ))

    Bnt=(0.035**2)    !  Use the Bnt given by Banner et al 2007
! ----------------------------------------------------------------------------------------------------------------
! #1        a1 and a2              get a1 and a2 as a function of U/cp
! ------------------------------------------------------------------------------------------------------------------

! find the fp and determine the cp

!          imax=maxloc(EDENS,1)
!          fpcheck=f(imax)
!          imax=-999

!          fp=fpinput

! needed: A(f), see Young and Babanin eq 19.
! Normally, it is read in and used, but per recommendation by Alex that it is unnecessary/obsolete, I set it to 1.0 here:
    ANAR=1.0 ! option 1
!   ANAR=ANAR_IN ! option 2

! (point output write location 1)
    if(.false.)then
       write(*,*)'a1,a2 = ',a1,a2
    endif
! two matrix operations follow:
    EDENST= (2.*g9**2) * ((2.0*pi)**(-4)) * (f**(-5)) * (ANAR**(-1)) * Bnt ! see eq 9 of "Implementation of new..."  Bnt is sigma_thr
    DEDENS=EDENS-EDENST
! if DEDENS < 0 set to zero

    DEDENS=MAX(0.0,DEDENS)


!    DO IS = 1, nfreq
!      WRITE(*,*) IS, EDENS(IS), EDENST(IS), DEDENS(IS)
!    END DO


!   ELIM=maxval(Edens)*1.0e-5 ! needed for v42D (can comment but don't delete)

    do is=1,nfreq

!v42U normalize by EDENST, a la Fabrice, this means that, with L=2 and large exceedence, the dissipation will get very strong
!v42U     NDEDENS(is)=DEDENS(is)/EDENST(is)
!v42D: normalize by EDENS, a la Kakha, this means that the normalized value is always less than 1, which means that using L=2 instead of L=1 will actually make the dissipation weaker

!             if(Edens(is).gt.ELIM)then
!                NDEDENS(is)=DEDENS(is)/EDENS(is) ! v42D
!             else
!                NDEDENS(is)=0.0
!             end if

       NDEDENS(is)=DEDENS(is)/EDENST(is) ! v42U

    enddo

    NDEDENS=MAX(0.0,NDEDENS)

! -------------------------------------------------------------------------------
!                      calculate T1
!  -------------------------------------------------------------------------------

!   T1=a1*f*A*DEDENS ! original matrix operation
    do  is=1,nfreq
       T1(is)=a1*f(is)*ANAR(is)*NDEDENS(is)**LPOW
    end do

! ------------------------------------------------------------------------------------
!                      calculate  T2  an integration from fp to f
! ------------------------------------------------------------------------------------

! OLD version: uses fp that falls on a node
! New version: provided peak does not necessarily fall on a node
! New New version: don't worry about fp, since stuff below fp is not breaking, so it isn't part of the calculation anyway.

    xff=0.
    do  is=1,nfreq
       ASUM=0.
       do ii=1,is
          xff(ii)=f(ii)
          ADF(ii)=ANAR(ii)*NDEDENS(ii)**MPOW
       enddo
       CALL integrate(ASUM,xff,ADF,is)
       T2(is)=a2*ASUM
    enddo

    T1=max(0.,T1)  ! WER a matrix operation
    T2=max(0.,T2)  ! WER a matrix operation
    Kds=(T1+T2)     ! WER a matrix operation

! (point output write location 2)

    if(.false.)then
       do  is=1,nfreq
          write(*,205)f(is),EDENS(is),EDENST(is),T1(is),T2(is)
       end do
    endif
205 format('threshold: ',5(1x,D12.5))

! figure out integrated T1 and T2 (for output purposes only)
! 3/9/2009 Since we have already non-dimensionalized by using NDEDENS instead of DEDENS, we do it backwards: Sds=Kds*Edens

    do is=1,nfreq
       Sds(is)=Kds(is)*Edens(is) ! note that in most cases, the calling routine will not use Sds
       ST1(is)=T1(is)*Edens(is)
       ST2(is)=T2(is)*Edens(is)
    end do

    DO IS = 1, MSC
      DO ID = 1, MDC
        SSDS(IS,ID) = KDS(IS)
        IF (ICOMP .GE. 2) THEN
          IMATDA(IS,ID) = IMATDA(IS,ID) + Kds(IS)
        ELSE
          IMATDA(IS,ID) = IMATDA(IS,ID) - kds(IS)
          IMATRA(IS,ID) = IMATRA(IS,ID) - kds(IS) * ACLOC(IS,ID)
        END IF
      END DO
    END DO

    CALL integrate(ST1_INT,f,ST1,nfreq)
    CALL integrate(ST2_INT,f,ST2,nfreq)
    CALL integrate(Sds_int,f,Sds,nfreq)

! (point output write location 3)
    if(.false.)then
       WRITE(*,*)'integral of T1,T2,Sds = ',ST1_INT,ST2_INT,Sds_INT
    endif

  END subroutine calc_Sds

! ************************************************************************************************

! -------------------------------------------------------------------------------------------------------------------------------
!       ====================================================

  SUBROUTINE integrate(ansNum, x,y,np)
!       ====================================================
    USE DATAPOOL, ONLY : RHOA, RHOW
    IMPLICIT NONE

!
!    Use Trapezoidal Rule to Integrate the area under the curve
!    specified by the data points in x and y
!
!    John Mahaffy 4/9/95
!
!
!   Local Variables
!
!    esterr   -   estimated error for Trapizoidal integration
!    sum2     -   Trapizoidal integration using every other available point
!

    REAL             , INTENT(IN)  ::  X(:)
    REAL             , INTENT(IN)  ::  Y(:)

    integer i
!  real sum2, esterr
    REAL          , INTENT(out) :: ANSNUM

    INTEGER           , INTENT(IN) :: NP

    ansnum=0.0
!  np=size(x) ! we cannot do this because array may be partially filled

    do  i=1,np-1
       ansnum=ansnum+.5*(y(i)+y(i+1))*(x(i+1)-x(i))
    end do

!
!   The following only works if np is an odd integer
!
!  sum2=0.0
!  do  i=1,np-1,2
!     sum2=sum2+.5*(y(i)+y(i+2))*(x(i+2)-x(i))
!  end do


    return
  end  SUBROUTINE integrate

subroutine calc_Lfactor(ip,Lfactor_S,S_in,DDIR_RAD,SIGMA_S,FRINTF,CINV_S,grav,WIND10,TESTFL)

    use datapool, only : rhoa, rhow, pi, ufric
    implicit none

! objective : provide Lfactor (1-d array) by solving for optimal REDUC (scalar). It uses a solver method
!     that should work for all monotonic relations. For non-monotonic relations, it may get confused.
!     Here, tau_normal definitely has a monotonic dependence on REDUC, so this potential limitation on
!     the solver is OK. (Monoticity of S_in, or lack thereof, is irrelevant.)
! points made here:
!  1) we have no idea what the tangential stress is, but we no that it isn't less
!              that zero, so we say that the normal stress cannot be greater than the total stress.
!  2) we don't know what the contribution is beyone fmax, but we know that it isn't negative, so again,
!              we  we say that the normal stress in our prognostic range cannot be greater than the
!              total stress...so tau_total=tau_normal_prognostic+tau_normal_diagnostic+tau_tangential...all are >=0,
!              so maximum allowed value of tau_normal_prognostic is tau_total
!  3) similar to Kakha's argument, we use the exponential adjustment that is stronger at the higher
!              frequencies (more reduction) since this is where our formula is less well-informed by
!              observations (more wiggle room), since typical observations are for the dominant waves, f=fp
!  4) Kakha uses f/fp to scale the reduction, which makes sense in the context of item (3), but use
!              of fp has disadvantages, so we use U/C instead. This means that for young wind sea,
!              S_in for the entire spectrum is reduced...will this be a problem? (the original thinking
!              was to use an fp value in a manner similar to Tolman and Chalikov, which is based on U/C somehow.
!  5) S_in=S_in_linear+S_in_exponential....same argument as (1)

    Integer, Intent(In)                           ::  ip
    REAL             , INTENT(OUT) ::  Lfactor_S(:)
    REAL             , INTENT(IN)  ::  S_in(:,:)
    REAL             , INTENT(IN)  ::  DDIR_RAD
    REAL             , INTENT(IN)  ::  FRINTF
    REAL             , INTENT(IN)  ::  CINV_S(:)
    REAL             , INTENT(IN)  ::  SIGMA_S(:)
    REAL             , INTENT(IN)  ::  grav
    REAL             , INTENT(IN)  ::  WIND10
    LOGICAL          , INTENT(INOUT)  ::  TESTFL

    REAL             , ALLOCATABLE  ::  S_in1D_S(:)
    REAL             , ALLOCATABLE  ::  S_in1D_L(:)
    REAL             , ALLOCATABLE  ::  DF(:)
    REAL             , ALLOCATABLE  ::  CINV_L(:)
    REAL             , ALLOCATABLE  ::  SIGMA_L(:)
    REAL             , ALLOCATABLE  ::  LFACTOR_L(:)

    REAL :: DSIGMA
    REAL :: REDUC
    REAL :: tau_normal
    REAL :: sign_new,sign_old
    REAL :: frat
    REAL :: err
    REAL :: RCHANGE
    REAL :: tau_total
    REAL :: TAUWLIM
    REAL :: U10MOD
    REAL :: freq_tmp
    REAL :: Cv,tauv

    INTEGER :: nf_old,nf_new
    INTEGER :: iter
    INTEGER :: slow_down
    INTEGER :: isol
    INTEGER :: IS

    TESTFL = .false.

    nf_old=size(S_in,1)

! Initialize LFACTOR_S (for error catching)
    LFACTOR_S = HUGE(LFACTOR_S)

    ALLOCATE(S_in1D_S(nf_old))

! To simplify calcs, use S_in1D(f)
    S_in1D_S = sum(S_in,1)*DDIR_rad ! note that variable is S_in(ID,IS)
    S_in1D_S = S_in1D_S*(2.0*PI) ! was in units of m2/(rad-Hz), so put in units m2/Hz

! find nf_new
    frat=FRINTF+1.0 ! =freq(nf_old)/freq(nf_old-1)
    FREQ_TMP = SIGMA_S(nf_old)/(2.0*pi)
    nf_new = -99
    do IS=(nf_old+1),(nf_old+100)
       FREQ_TMP=FREQ_TMP*frat
       if(FREQ_TMP > 10.0)then ! 10.0 here is f=10 Hz
          nf_new=IS
          exit
       endif
    enddo

! allocate arrays on nf_new
    ALLOCATE(S_in1D_L(nf_new))
    ALLOCATE(DF(nf_new)) ! REQ'D?
    ALLOCATE(CINV_L(nf_new))
    ALLOCATE(SIGMA_L(nf_new))
    ALLOCATE(LFACTOR_L(nf_new))

! extend freqs to nf_new
    SIGMA_L=SIGMA_S
    CINV_L=CINV_S
    S_in1D_L=S_in1D_S
    do IS=(nf_old+1),nf_new
       SIGMA_L(IS)=SIGMA_L(IS-1)*frat
       CINV_L(IS)=SIGMA_L(IS)*0.102 ! deep water assumption for hf tail
! For Sin, we use an f-2 approximation...Sin=Beta*E~sigma*(U/C)^2*E=f1*f2*f-5=f-2
       S_in1D_L(IS)=S_in1D_L(nf_old)*(SIGMA_L(nf_old)/SIGMA_L(IS))**2
    enddo

    tau_total=(UFRIC(IP)**2)*rhoa
    U10MOD=UFRIC(IP)*28.0
    Cv=(-5.0e-5)*WIND10+1.1e-3
    Cv=max(Cv,0.0) ! otherwise Cv<0 for U10>22 m/s
    tauv=rhoa*Cv*(WIND10**2)
    TAUWLIM=tau_total-tauv

    DO IS=1,NF_new
       DSIGMA=SIGMA_L(IS)*FRINTF
       DF(IS)=DSIGMA/(2.0*PI)
    END DO

! calculate without reduction
    REDUC=0.0
    call calc_tau_normal(tau_normal,Lfactor_L,REDUC,S_in1D_L,df,CINV_L,U10MOD,grav,rhow)

    if(TESTFL)then
206    format('tau_total = ',f9.6,' ; tauv = ',f9.6,' ; tau_normal = ',f9.6,' ; TAUWLIM = ',f9.6)
       write(*,206)tau_total,tauv,tau_normal,TAUWLIM
    endif

    if(tau_normal < TAUWLIM)then
       DO IS=1,NF_OLD
          Lfactor_S(IS)=1.0
       END DO
       RETURN
    else
       if(TESTFL)then
          write(*,*)'reducing Sin to make tau_normal match TAUWLIM'
       endif
    endif

!    if(tau_normal > 10.0) then
!      write(*,*) ip, wind10, ufric, tau_normal 
!    end if

! calculate with reduction, start value arbitrary

    REDUC=0.05
    call calc_tau_normal(tau_normal,Lfactor_L,REDUC,S_in1D_L,df,CINV_L,U10MOD,grav,rhow)
    err=tau_normal-TAUWLIM
    sign_new=sign(1.0,err)
    slow_down=0
    RCHANGE=2.0
    isol=0

    do iter=1,500

       if(tau_normal > TAUWLIM)then ! increase REDUC
          REDUC=REDUC*RCHANGE
       else ! then reduce REDUC
          REDUC=REDUC/RCHANGE
       endif

       call calc_tau_normal(tau_normal,Lfactor_L,REDUC,S_in1D_L,df,CINV_L,U10MOD,grav,rhow)

       err=tau_normal-TAUWLIM

       sign_old=sign_new

       sign_new=sign(1.0,err)

       if(TESTFL)then
207       format('iter = ',i3,' ; REDUC = ',f9.6,' ; RCHANGE = ',f9.6,' ; tau_normal = ',f9.6,' ; err = ',f9.6)
          write(*,207)iter,REDUC,RCHANGE,tau_normal,err
       endif

       if(abs(sign_new-sign_old) .gt. tiny(1.))then
         RCHANGE=0.5*(1.0+RCHANGE)
       endif

       if((abs(err)/TAUWLIM) < 1.656e-06)then
          isol=1
          exit
       endif

    end do

    if(isol==0)then
       write(*,*)'no solution found'
       write(*,*)'  tau_normal = ',tau_normal,' TAUWLIM = ',TAUWLIM,' err = ',err,'(abs(err)/TAUWLIM)  = ',(abs(err)/TAUWLIM)
       write(*,*) wind10, ufric
       if((abs(err)/TAUWLIM).ge.1.0)then
          open(401,file='debug.dat',form='unformatted')
          write(401)nf_new     ! scalar
          write(401)Lfactor_L  ! LFACTOR_L(nf_new)
          write(401)S_in1D_L   ! S_in1D_L(nf_new)
          write(401)df         ! DF(nf_new)
          write(401)CINV_L     ! CINV_L(nf_new)
          write(401)TAUWLIM    ! scalar
          write(401)tau_normal ! scalar
          write(401)REDUC      ! scalar
          write(401)U10MOD     ! scalar
          write(401)rhow       ! scalar
          close(401)
          stop
       end if
    end if

    DO IS=1,NF_OLD
       Lfactor_S(IS)=Lfactor_L(IS)
    END DO

  end subroutine calc_Lfactor

  subroutine calc_tau_normal(tau_normal,Lfactor,REDUC,S_in1D,df,CINV,U10MOD,grav,rhow)
! objective : apply exponential reduction.. more reduction at higher U/C with factor smoothly intersecting value=1 at U/C=1. It is not applied for U/C<1 since that would produce an increase...so factor 1 for U/C<=1

    implicit none

    REAL             , INTENT(OUT) ::  tau_normal
    REAL             , INTENT(OUT) ::  Lfactor(:)
    REAL             , INTENT(IN)  ::  REDUC
    REAL             , INTENT(IN)  ::  S_in1D(:)
    REAL             , INTENT(IN)  ::  df(:)
    REAL             , INTENT(IN)  ::  CINV(:)
    REAL             , INTENT(IN)  ::  U10MOD
    REAL             , INTENT(IN)  ::  grav
    REAL             , INTENT(IN)  ::  rhow

    REAL             , ALLOCATABLE  :: A(:)
    REAL             , ALLOCATABLE  :: UoverC(:)
    REAL             , ALLOCATABLE  :: S_in1D_red(:)

    INTEGER :: nf
    INTEGER :: IS

    nf=size(df)

    ALLOCATE(UoverC(nf))
    ALLOCATE(A(nf))
    ALLOCATE(S_in1D_red(nf))

    UoverC=U10MOD * CINV ! matrix op

    A=exp((1-UoverC) * REDUC) ! matrix op
    do IS=1,nf
       Lfactor(IS)=min(1.0,A(IS)) ! cannot be done as matrix op?
    enddo
    S_in1D_red=S_in1D * Lfactor ! matrix op

    tau_normal=0.0
    do IS=1,nf
       tau_normal=tau_normal+rhow*grav*S_in1D_red(IS)*CINV(IS)*df(IS)
    enddo

  end subroutine calc_tau_normal


!****************************************************************
  SUBROUTINE SSWELL (IP,ETOT,ACLOC,IMATRA, IMATDA,URMSTOP,CGo)

! note that CGo is for diagnostic purposes only, may be excluded
! excluded : SPCDIR AC2 DEP2 IMATRA

!****************************************************************
    USE DATAPOOL
    IMPLICIT NONE


    INTEGER, INTENT(IN)    :: IP
    REAL   , INTENT(IN)    :: ETOT
    REAL   , INTENT(IN)    :: CGo(:,:)   ! CGo(MSC,MICMAX) ! group velocity without currents, but includes depth effects
    REAL   , INTENT(IN)    :: URMSTOP

    REAL   , INTENT(INOUT) :: IMATDA(:,:) ! IMATDA(MDC,MSC)
    REAL   , INTENT(INOUT) :: IMATRA(:,:) ! IMATDA(MDC,MSC)
    REAL   , INTENT(IN)    :: ACLOC(:,:) ! ACLOC(MDC,MSC)

!   local constants:
    REAL, parameter   :: feswell=0.0035 ! "fe is of the order 0.002 to 0.008" (Fabrice proposes 0.0035)
    REAL, parameter   :: Re_crit=1.0e+5 ! Re_c or SWELLF4
    REAL, parameter   :: Cdsv=1.2 ! C_dsv or SWELLF5
    REAL, parameter   :: nu_air=1.4E-5 

!   local variables:
    INTEGER           :: IS
    INTEGER           :: ID
    REAL              :: aorb,USIGTOP
    REAL              :: Re
    REAL              :: SWDIS(MSC) ! this is S_SWELL / E(f,theta)
!    REAL              :: myu,SWDIS_CHECK ! error checking, may be omitted

! Method:
! In general....
! @N/@t=S/sig=C(sig)*E(sig,theta)/sig=C(sig)*N(sig,theta)
! where C(sig) has units of rad/sec
! ...so if we are passing to IMATDA, we just need to pass C(sig)
! Here, SWDIS(IS) is my C(sig)

    aorb=2.0*sqrt(ETOT)
    USIGTOP=URMSTOP*SQRT(2.0)
    Re=4.0*USIGTOP*aorb/nu_air
    IF(Re > Re_crit)then
       DO IS=1, MSC 
          SWDIS(IS) = RHOAW * 16.0 * feswell * SPSIG(IS)**2  & ! note that "-1" omitted since LHS
          &           * URMSTOP / g9 ! units of radian^2/s (radians/sec is conventional)
       END DO
    else
       DO IS=1, MSC 
          SWDIS(IS) = RHOAW * Cdsv * 2.0 * WK(IP,IS)  & ! note that "-1" omitted since LHS
          &           * SQRT(2.0 * NU_AIR * SPSIG(IS)) ! units of radian^1.5/s (radians/sec is conventional)

! start error check block (may be omitted)
!                myu=2.0*RHOAW*(1.0/g9)*(1.0/CGo(IS,1))
!     &              *(SPCSIG(IS)**2.5)*sqrt(2.0*NU_AIR)
!                write(*,*)'myu(',IS,') = ',myu
!                SWDIS_CHECK=Cdsv*myu*CGo(IS,1)
!                write(*,*)'SWDIS(IS),SWDIS_CHECK = ',
!     &                    SWDIS(IS),SWDIS_CHECK
! end error check block (may be omitted)

       END DO
    endif

    DO IS=1, MSC 
      DO ID = 1, MSC
         IMATDA(IS,ID) = IMATDA(IS,ID) + SWDIS(IS)
         IF (ACLOC(IS,ID) .GT. THR) IMATRA(IS,ID) = SWDIS(IS) * ACLOC(IS,ID)
      END DO
    END DO

    RETURN
  END SUBROUTINE SSWELL

  SUBROUTINE SWIND_DBYB (IP,WIND10,THETAW,IMATRA,SSINE)   ! 30.21

    use datapool
    implicit none

    INTEGER, INTENT(IN) :: IP

 
    INTEGER   :: ID, IS

    REAL      :: THETAW, SIGMA, SWINEB, CTW, STW, COSDIF, TEMP2, UoverC, WIND10
!
    REAL      :: S_IN(MSC,MDC), IMATRA(MSC,MDC)
    REAL, INTENT(OUT) :: SSINE(MSC,MDC)
!
    REAL      :: DMAX,AINV
    REAL      :: KTHETA(MSC,MDC) ! like D(theta), except max value at each freq is unity
    REAL      :: ASPR(MSC), SIGDENS(MSC), SQRTBN(MSC), CINV(MSC), LFACTOR(MSC)
    REAL      :: GAMMAD, GDONEL, TEMP4, TEMP42, TEMP5, TEMP6, BN
    LOGICAL   :: TESTFL

!WER      REAL HM0,ETOT,EDENS2D ! FOR test calcs (remove later)

    CTW   = COS(THETAW)                                       !          40.41
    STW   = SIN(THETAW)                                       !          40.41

    IF (WIND10 .LT. VERYSMALL) RETURN

    DO  IS = 1, MSC
       SIGDENS(IS) = 0.
       DO  ID = 1, MDC
          SIGDENS(IS) = SIGDENS(IS) + SPSIG(IS) * AC2(IP,IS,ID)
          KTHETA(IS,ID) = AC2(IP,IS,ID)
       END DO
       SIGDENS(IS)=SIGDENS(IS)*DDIR 
    END DO

!    DO IS = 1, MSC
!       DMAX=-1.0
!       DO ID = 1, MDC
!         IF(KTHETA(IS,ID).GT.DMAX)DMAX=KTHETA(IS,ID)
!       END DO
!       IF(DMAX.EQ.0.0)THEN ! FIX FOR FREQ BINS (USUALLY FIRST TWO OR SO) THAT ARE EMPTY
!          DMAX=1.0
!          DO  ID = 1, MDC
!             KTHETA(IS,ID)=1.0
!          END DO
!       ENDIF
!       DO ID = 1, MDC
!         KTHETA(IS,ID)=KTHETA(IS,ID)/DMAX
!       END DO
!    END DO

    DO IS = 1, MSC
      DMAX = 0.
      DMAX = MAXVAL(KTHETA(IS,:))
      IF (DMAX .LT. 10.E-10) THEN
        KTHETA(IS,:) = 1.
      ELSE
        KTHETA(IS,:)=KTHETA(IS,:)/DMAX
      END IF
    END DO

    DO IS = 1, MSC
       AINV=0.0
       DO ID = 1, MDC
          AINV=AINV+KTHETA(IS,ID)*DDIR
       END DO
       ASPR(IS)=1./AINV
       BN=ASPR(IS)*SPSIG(IS)**5*SIGDENS(IS)/(2.*G9**2)
       SQRTBN(IS)=SQRT(BN)
    END DO

    TEMP2 = 28.0 * UFRIC(IP)
    S_IN=0.0
    DO IS = 1, MSC
       SIGMA = SPSIG(IS)                                  
       CINV(IS)  = WK(IP,IS) / SIGMA
       UoverC = TEMP2 * CINV(IS)
       DO ID=1,MDC
          IF ( WIND10 .GT. THR ) THEN
             COSDIF = COSTH(ID)*CTW+SINTH(ID)*STW          
             TEMP4=( UoverC * COSDIF - 1.0 )
             TEMP4=MAX(0.0,TEMP4)
             TEMP42=TEMP4**2
             TEMP5=10.0*SQRTBN(IS)*TEMP42-11.0
             TEMP6=1.0+TANH(TEMP5)
             GDONEL=2.8-TEMP6
             GAMMAD=GDONEL*SQRTBN(IS)*TEMP42
             SWINEB =  GAMMAD * SIGMA * RHOAW ! NEW SWINEB
             S_IN(IS,ID)=SWINEB*AC2(IP,IS,ID)*SPSIG(IS)
          END IF
       ENDDO
    ENDDO

    call calc_Lfactor(ip,Lfactor,S_in,DDIR,SPSIG,FRINTF,CINV,G9,WIND10,TESTFL)

    SSINE = S_IN

    DO IS = 1, MSC 
       DO ID = 1, MDC
          IF ( WIND10 .GT. VERYSMALL ) THEN
             S_IN(IS,ID)= S_IN(IS,ID) * Lfactor(IS)
             IMATRA(IS,ID) = IMATRA(IS,ID) + S_IN(IS,ID)/SPSIG(IS) 
             SWINEB=S_IN(IS,ID)/(AC2(IP,IS,ID)*SPSIG(IS))
          END IF
       ENDDO
    ENDDO

  END SUBROUTINE SWIND_DBYB

END MODULE SdsBabanin
