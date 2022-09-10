#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This can remain ...
      SUBROUTINE PEAK_PARAMETER_BREAK(IP,WALOC,LPP)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)           :: IP
         REAL(rkind), INTENT(IN)       :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)      :: LPP
         REAL(rkind)                   :: KPP, FPP, CPP, WNPP, CGPP, TPP

         INTEGER                       :: IS, ID, ISIGMP
         REAL(rkind)                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, WKDEPD, WNPD
         REAL(rkind)                   :: DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, CPWN
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(WALOC)

       IF (MAXAC .gt. VERYSMALL .AND.  DEP(IP) .GT. DMIN) THEN

         ETOTF3 = ZERO
         ETOTF4 = ZERO
         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             HQUOT  = WALOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
           END DO
         END DO

         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN
           FPP    = ETOTF3/ETOTF4*PI2
           CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           TPP    = PI2/FPP
           LPP    = PI2/KPP

         ELSE 
           LPP = ZERO
         END IF
       ELSE
         LPP = ZERO
       END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, WALOC, SSBR, DSSBR)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: IP

      REAL(rkind), INTENT(IN)   :: WALOC(NUMSIG,NUMDIR), SME, KME, ETOT, HS

      REAL(rkind), INTENT(INOUT)     :: SSBR(NUMSIG,NUMDIR), DSSBR(NUMSIG,NUMDIR)
      REAL(rkind)   :: LPP

      REAL(rkind) :: BETA, QQ, QB, BETA2, ARG
      REAL(rkind) :: S0
#ifdef DEBUG
      integer, save :: idxcall = 0
#endif  
      REAL(rkind) :: AUX, GAMMA_WB, COEFF_A 
      REAL(rkind) :: SBRD, WS, SURFA0, SURFA1, COEFF_B, SURFSEL

      INTEGER :: IS, ID
      REAL(rkind) :: BJALFA, SDBC1, CBJ, TM01, TM02, FMEANloc
	  
!     Compute breaking fraction !
      SELECT CASE(ICRIT)
       CASE(1) ! simple breaking coefficient
         HMAX(IP) = BRHD * DEP(IP)
       CASE(2) ! Suggestion of Dingemans 
         IF (KME .GT. VERYSMALL) THEN
           S0    = HS / (PI2/KME) 
           GAMMA_WB  = 0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
           HMAX(IP)  = GAMMA_WB * DEP(IP)
         ELSE
           HMAX(IP)  = BRHD * DEP(IP)
         END IF
       CASE(3) ! D. based on peak steepness 
         CALL PEAK_PARAMETER_BREAK(IP,WALOC,LPP)
         IF (LPP .GT. VERYSMALL) THEN
           S0    = HS/LPP
           GAMMA_WB =  0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
           HMAX(IP) = GAMMA_WB * DEP(IP)
         ELSE
           HMAX(IP) = BRHD * DEP(IP)
         END IF
       CASE DEFAULT
         CALL WWM_ABORT('ICRIT HAS A WRONG VALUE')
      END SELECT
	  
!     Variable breaker index
      CALL COMPUTE_VAR_GAMMA(IP)
	  
!     Transform to monochromatic waves !
      IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(TWO)
	  
!     Compute beta ratio ! 
      IF ( (HMAX(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
		IF (IBREAK == 1) THEN
          BETA = SQRT(8.d0*ETOT)/HMAX(IP)
        ELSEIF (IBREAK == 3) THEN
          !BETA = SQRT(8.d0*ETOT)/(BRHD/sqrt(2.d0)*DEP(IP)) 
          BETA = SQRT(8.d0*ETOT)/(GAMMA(IP)*DEP(IP)) 
        ENDIF
        BETA2 = BETA**2.d0
      ELSE
        BETA = ZERO 
        BETA2 = ZERO 
      END IF

      IF (BETA <= 0.5_RKIND) THEN
        QQ = ZERO 
      ELSE IF (BETA <= ONE) THEN
        QQ = (TWO*BETA-ONE)**2
      END IF
	  
!     Compute breaking fraction based on the idea of Henrique Alves !
      IF ( BETA .LT. 0.2_rkind ) THEN
        QB     = ZERO
      ELSE IF ( BETA .LT. ONE ) THEN
        ARG    = EXP  (( QQ - ONE ) / BETA2 )
        QB     = QQ - BETA2 * ( QQ - ARG ) / ( BETA2 - ARG )
        DO IS = 1, 3
          QB     = EXP((QB-ONE)/BETA2)
        END DO
      ELSE
        QB = ONE - SMALL
      ENDIF
! 
      QBLOCAL(IP) = QB
!
      IF (IBREAK == 1) THEN ! Battjes & Janssen
        IF ( BETA2 .GT. SMALL  .AND. MyABS(BETA2 - QB) .GT. SMALL ) THEN
          IF ( BETA2 .LT. ONE - SMALL) THEN
            SURFA0  = - ( ALPBJ / PI) *  QB * SME / BETA2
          ELSE
            SURFA0  = - (ALPBJ/PI) * SME
          END IF
        ELSE
          SURFA0 = ZERO
        END IF
        SURFA1 = SURFA0 
      ELSEIF (IBREAK == 2) THEN ! Battjes & Janssen SWAN code works only for implicit since it is linearized and both terms are positive the explanation is given in the SWAN code  
        IF ( BETA2 .GT. SMALL  .AND. MyABS(BETA2 - QB) .GT. SMALL ) THEN
          IF ( BETA2 .LT. ONE - SMALL) THEN
            WS   = (ALPBJ / PI) *  QB * SME / BETA2
            SbrD = WS * (ONE - QB) / (BETA2 - QB)
          ELSE
            WS   = (ALPBJ/PI) * SME
            SbrD = ZERO 
          END IF
          SURFA0 = SbrD ! right hand side is positive to balance left hand side underelax by SbrD ! explicit schemes cannot works since this term is positive!!!
          SURFA1 = - WS + SbrD ! will be inverted in the implicit integration ... 
        ELSE
          SURFA0 = ZERO 
          SURFA1 = ZERO 
        END IF
      ELSEIF (IBREAK == 3) THEN ! Thornton & Guza 1983
        IF (ICOMP .GE. 2) THEN ! Implicit
          IF ( BETA2 .GT. 0D0 ) THEN
!            COEFF_A = 0.42_rkind ! TG83 default (Torrey Pine beach)
            COEFF_A = BRHD/sqrt(2.d0) !Guerin
            COEFF_B = 4.0_rkind
            IF ( BETA2 .LT.1D0 ) THEN
              WS   = 75D-2*COEFF_A*ALPBJ**3.d0*SME*BETA2**(0.5d0*(COEFF_B+1.0_rkind))/MyREAL(SQRT(PI))
              SbrD = 5D-1*MyREAL(3.d0+COEFF_B)*WS
            ELSE
              WS   = 75D-2*COEFF_A*ALPBJ**3.d0*SME/MyREAL(SQRT(PI))
              SbrD = WS
            ENDIF
            SURFA0 = SbrD - WS
            SURFA1 = SbrD
          ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
          ENDIF 
        ELSE ! Explicit
          IF ( BETA2 .GT. 0D0) THEN
		    !QB = (SQRT(8.d0*ETOT)/(BRHD/sqrt(2.d0)*DEP(IP)))**4 
		    QB = (SQRT(8.d0*ETOT)/(GAMMA(IP)*DEP(IP))) 
			IF ( QB .LT. 1D0 ) THEN
              SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*QB**4*SQRT(8.d0*ETOT)/DEP(IP)
			ELSE
              SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
			ENDIF
          ELSE
            SURFA0 = 0D0
          ENDIF
        ENDIF
      ELSEIF (IBREAK == 4) THEN ! Church and Thornton (1993)
        IF (ICOMP .LT. 2) THEN ! Only explicit for now
          QB = (SQRT(8.d0*ETOT)/(GAMMA(IP)*DEP(IP)))
		  SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
			     & *(1+tanh(8*(QB-1)))*(1-(1+QB**2)**(-5/2))
		  !SURFA0 = -(0.75d0*BTG(IP)**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
			     !& *(1+tanh(8*(QB-1)))*(1-(1+QB**2)**(-5/2))
        ENDIF
      ELSEIF (IBREAK == 5) THEN ! Baldock et al. (1998)
        IF ( BETA2 .GT. 0D0) THEN
		  QB = exp(-(GAMMA(IP)*DEP(IP))/(SQRT(8.d0*ETOT))) 
		  IF ( QB .LT. 1D0 ) THEN
            SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*QB**4*SQRT(8.d0*ETOT)/DEP(IP)
		  ELSE
            SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
		  ENDIF
        ELSE
          SURFA0 = 0D0
        ENDIF
        !IF (ICOMP .LT. 2) THEN ! Only explicit for now
		!  IF(GAMMA(IP).LE.0.2) GAMMA(IP) = 0.40
        !  QB = (GAMMA(IP)*DEP(IP))/(SQRT(8.d0*ETOT))
		!  SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
		!	     & *(1 + (4/(3*MyREAL(PI)))*(QB**3+1.5D0*QB)*exp(-QB**2) - erf(Qb))
        !ENDIF
      ENDIF
	  
!     Copy Right hand side and diagonal term 
      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          DSSBR(IS,ID)  = SURFA1
          SSBR(IS,ID)   = SURFA0 * WALOC(IS,ID)
        END DO
      END DO 

#ifdef SCHISM
      SBR(:,IP) = ZERO
      DO IS=1,NUMSIG
        DO ID=1,NUMDIR
          SBR(1,IP)=SBR(1,IP)+G9*COSTH(ID)*WK(IS,IP)*SSBR(IS,ID)*DS_INCR(IS)*DDIR !TG: SPSIG(IS) is simplified
          SBR(2,IP)=SBR(2,IP)+G9*SINTH(ID)*WK(IS,IP)*SSBR(IS,ID)*DS_INCR(IS)*DDIR !TG: SPSIG(IS) is simplified
        ENDDO
      ENDDO
#endif
      END SUBROUTINE 
!**********************************************************************
!*  Computing the local breaker index depending on the local slope and wave info
!**********************************************************************
      SUBROUTINE COMPUTE_VAR_GAMMA(IP)
       USE DATAPOOL
	   USE schism_glbl, only : dldxy,indel
       IMPLICIT NONE

       INTEGER, INTENT(IN) :: IP
	   INTEGER :: ICOUNT, IE, INNE, INDCUT
	   
       REAL(rkind) :: HS, TM01, TM02, TM10, KLM, WLM, ETOTS, ETOTC, DM, DSPR, TMPDIR
       REAL(rkind) :: FPP, TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKDSPR, PEAKD, DPEAK, TPPD, KPPD, CGPD, CPPD
       REAL(rkind) :: GAMMA_BETA, GAMMA_KD, GAMMA_0, a1, a2, a3
	   
	   ! Computing the mean wave characteristics
       IF (OUTFRHIGHFLAG) THEN
         INDCUT = MINLOC(ABS(SPSIG/PI2-OUTFRHIGH),1)
         IF (SPSIG(INDCUT)/PI2.LT.OUTFRHIGH) THEN
           INDCUT = INDCUT + 1
         ENDIF
         !CALL MEAN_PARAMETER(IP,AC2(:,:,IP),INDCUT,HS,TM01,TM02,TM10,KLM,WLM)
		 CALL MEAN_DIRECTION_AND_SPREAD(IP,AC2(:,:,IP),INDCUT,ETOTS,ETOTC,DM,DSPR)
         CALL PEAK_PARAMETER(IP,AC2(:,:,IP),INDCUT,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)
       ELSE
         !CALL MEAN_PARAMETER(IP,AC2(:,:,IP),NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
		 CALL MEAN_DIRECTION_AND_SPREAD(IP,AC2(:,:,IP),NUMSIG,ETOTS,ETOTC,DM,DSPR)
         CALL PEAK_PARAMETER(IP,AC2(:,:,IP),NUMSIG,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)
       ENDIF
	   
	   ! Computation of the bed slope in the mean wave direction
       TMPDIR = DM*PI/180.d0 + PI/2.d0 ! Conversion from nautical to math. and in rad.
       TANBETA_WD(IP) = TANBETA_X(IP)*cos(TMPDIR) + TANBETA_Y(IP)*sin(TMPDIR)
	   
	   ! Breaker index: Salmon et al. (2015)
	   GAMMA_0 = 0.35D0
	   a1 = 6D0  !7.59D0
	   a2 = 0D0  !-8.02D0
	   a3 = 5D0  !8.09D0
	   GAMMA_BETA = GAMMA_0 + a1*TANBETA_WD(IP)
	   GAMMA_KD   = a2 + a3*DEP(IP)*KLM
	   GAMMA(IP) = 0.50 + 4*TANBETA_X(IP)
	   !GAMMA(IP)  = GAMMA_BETA/tanh(GAMMA_BETA/GAMMA_KD)
	   !GAMMA(IP) = 0.65 + 4*TANBETA_WD(IP)/(DEP(IP)*KLM)
	   !GAMMA(IP) = 0.32 + 30*TANBETA_WD(IP)/(DEP(IP)*KPP)
	   BTG(IP)  = 1.2/(GAMMA(IP)*DEP(IP))+0.1
	   IF(DEP(IP).LE.2D0) BTG(IP) = 0.5
	   IF(BTG(IP).GE.1.5D0) BTG(IP) = 1.5D0
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************