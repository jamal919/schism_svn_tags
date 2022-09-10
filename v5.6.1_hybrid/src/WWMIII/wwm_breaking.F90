#include "wwm_functions.h"
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

      REAL(rkind) :: BETA, QB, BETA2, ARG
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
       CASE(4) ! Variable breaker index based on the beach slope
         CALL COMPUTE_VAR_GAMMA(IP)
         HMAX(IP) = GAMMA(IP) * DEP(IP)
       CASE DEFAULT
         CALL WWM_ABORT('ICRIT HAS A WRONG VALUE')
      END SELECT
	  
!     Transform to monochromatic waves !
      IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(TWO)
	  
	  IF (ICOMP .GE. 2) CALL WWM_ABORT('Use explicit schemes for breaking')
      SELECT CASE(IBREAK)
       CASE(1) ! Thornton & Guza 1983
		 QB = (SQRT(8.d0*ETOT)/HMAX(IP)) 
		 IF ( QB .LT. 1D0 ) THEN
           SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*QB**4*SQRT(8.d0*ETOT)/DEP(IP)
		 ELSE
           SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
		 ENDIF
       CASE(2) ! Church and Thornton (1993)
         BETA = SQRT(8.d0*ETOT)/HMAX(IP)
		 IF ( QB .LT. 1D0 ) THEN
           SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
		          & *(1+tanh(8*(BETA-1)))*(1-(1+BETA**2)**(-5/2))
		 ELSE
           SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
		 ENDIF
       CASE(3) ! Baldock et al. (1998)
         !QB = exp(-HMAX(IP)/SQRT(8.d0*ETOT)) 
		 !IF ( QB .LT. 1D0 ) THEN
         !  SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*QB**4*SQRT(8.d0*ETOT)/DEP(IP)
		 !ELSE
         !  SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
		 !ENDIF
		 BETA = HMAX(IP)/SQRT(8.d0*ETOT) ! creates a sigsev fault, who knows why...
         SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
			    & *(1 + (4/(3*MyREAL(sqrt(PI))))*(BETA**3+1.5*BETA)*exp(-BETA**2) - erf(BETA))
       CASE DEFAULT
         CALL WWM_ABORT('IBREAK HAS A WRONG VALUE')
      END SELECT
	  
!     Copy Right hand side and diagonal term 
      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          !DSSBR(IS,ID)  = SURFA1
          DSSBR(IS,ID)  = 0D0 ! only explicit now
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
       USE schism_glbl, only : tanbeta_x
       IMPLICIT NONE

       INTEGER, INTENT(IN) :: IP
	   INTEGER :: ICOUNT, IE, INNE, INDCUT
	   
       REAL(rkind) :: HS, TM01, TM02, TM10, KLM, WLM, ETOTS, ETOTC, DM, DSPR, TMPDIR
       REAL(rkind) :: FPP, TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKDSPR, PEAKD, DPEAK, TPPD, KPPD, CGPD, CPPD
       REAL(rkind) :: dpdx_e, dpdy_e, dpdx, dpdy, tmpx, tmpy
	   
	   ! Computing the mean wave characteristics
       !IF (OUTFRHIGHFLAG) THEN
         !INDCUT = MINLOC(ABS(SPSIG/PI2-OUTFRHIGH),1)
         !IF (SPSIG(INDCUT)/PI2.LT.OUTFRHIGH) THEN
           !INDCUT = INDCUT + 1
         !ENDIF
         !CALL MEAN_PARAMETER(IP,AC2(:,:,IP),INDCUT,HS,TM01,TM02,TM10,KLM,WLM)
		 !CALL MEAN_DIRECTION_AND_SPREAD(IP,AC2(:,:,IP),INDCUT,ETOTS,ETOTC,DM,DSPR)
         !CALL PEAK_PARAMETER(IP,AC2(:,:,IP),INDCUT,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)
       !ELSE
         !CALL MEAN_PARAMETER(IP,AC2(:,:,IP),NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
		 !CALL MEAN_DIRECTION_AND_SPREAD(IP,AC2(:,:,IP),NUMSIG,ETOTS,ETOTC,DM,DSPR)
         !CALL PEAK_PARAMETER(IP,AC2(:,:,IP),NUMSIG,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)
       !ENDIF
	   
	   ! Computation of the bed slope in the mean wave direction
       !TMPDIR = DM*PI/180.d0 + PI/2.d0 ! Conversion from nautical to math. and in rad.
       !TANBETA_WD(IP) = TANBETA_X(IP)*cos(TMPDIR) + TANBETA_Y(IP)*sin(TMPDIR) !Computation of the bed slope at node IP
	   
	   ! Breaker index
	   GAMMA(IP) = 0.55 + abs(2.8*tanbeta_x(IP))
	   !GAMMA(IP) = 0.42 + 2.5*tanbeta_x(IP)
	   !GAMMA(IP) = 0.32 + 30*TANBETA_WD(IP)/(DEP(IP)*KPP)
	   
	   ! B parameter (not used for now)
	   !BTG(IP)  = 1.2/(GAMMA(IP)*DEP(IP))+0.1
	   !IF(DEP(IP).LE.2D0) BTG(IP) = 0.5
	   !IF(BTG(IP).GE.1.5D0) BTG(IP) = 1.5D0
      END SUBROUTINE
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
	  REAL FUNCTION erf(x)
	     ! Error function from Numerical Recipes.
	     ! erf(x) = 1 - erfc(x)
         IMPLICIT NONE
         REAL :: dumerfc, x
         REAL :: t, z

         z = abs(x)
         t = 1.0 / ( 1.0 + 0.5 * z )

         dumerfc = t * exp(-z * z - 1.26551223 + t *	    &
	           ( 1.00002368 + t * ( 0.37409196 + t *		&
               ( 0.09678418 + t * (-0.18628806 + t *		&
			   ( 0.27886807 + t * (-1.13520398 + t *		&
               ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

         if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc     
	     erf = 1.0 - dumerfc
		 RETURN
	  END FUNCTION
