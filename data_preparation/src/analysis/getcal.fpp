C     =================
      SUBROUTINE GETCAL(vtx)
C     =================
C     This routine determines the total energy in the different parts
C     of the calorimeter, the transverse E_t in the event and the
C     maximum rapidity for a condensate with at least 400 MeV energy.
C     In the second part of this routine the calorimeter timing is
C     reconstructed.
C
C     A.Quadt & G.Howell, May 1997.
C     ----------------------------------------------------------------

      IMPLICIT NONE

#include "partap.inc"
#include "caltru.inc"
#include "cconsa.inc"
#include "const.inc"
c#include "ctrlcm.inc"
#include "zrevt.inc"
#include "ctime2.inc"
#include "sidat95.inc"
#include "common.inc"

      INTEGER   NO_USE, I, Ierr, vector(5), j
      DATA      NO_USE / -999.99 /
      REAL      RAD
      REAL      Z, TANGLE, ANGLE, ETA_CON
      REAL      HADE, HADPZ
      CHARACTER KIND*5      
      REAL      Ecut,ImBcut,Vtx(3),Lateness,Avgt(0:5), AvgTErr(0:5),
     +          TChiSq(0:5),TProb(0:5),ESum(0:5)
      INTEGER   NPMT(0:5)
      logical   usecell
      real thres
      data thres /0.4/
 
      real xv(3)
      real xpos, ypos, zpos, x1, y1, z1, r1, th1, ph1
      DATA      imbcut/0.35/
      SAVE      imbcut

      cal_tot = 0.0
      FEMC_EN  = 0.0
      BEMC_EN  = 0.0
      REMC_EN  = 0.0
      FHAC_EN  = 0.0
      BHAC_EN  = 0.0
      RHAC_EN  = 0.0
      ET_TOT = 0.0
      PZ_TOT2 = 0.0
      HADE = 0.0
      HADPZ = 0.0

      FCAL_TM  = NO_USE
      RCAL_TM  = NO_USE
      BCAL_TM  = NO_USE
      GLOB_TM  = NO_USE
      ETAMAX   = NO_USE


C --- CAL ENERGIES ---
C --------------------
      call setvtx(xv)
      DO I=1,COUTAB(CALTRU)

        CALL FETTAB(CALTRU, ID, I)
        cal_tot = cal_tot + caltru_E

        CALL CCWHAT (CALTRU_CELLNR, KIND, VECTOR, IERR)
        IF (KIND(1:2).EQ.'RE') REMC_EN = REMC_EN + CALTRU_E
        IF (KIND(1:2).EQ.'BE') BEMC_EN = BEMC_EN + CALTRU_E
        IF (KIND(1:2).EQ.'FE') FEMC_EN = FEMC_EN + CALTRU_E
        IF (KIND(1:2).EQ.'RH') RHAC_EN = RHAC_EN + CALTRU_E
        IF (KIND(1:2).EQ.'BH') BHAC_EN = BHAC_EN + CALTRU_E
        IF (KIND(1:2).EQ.'FH') FHAC_EN = FHAC_EN + CALTRU_E
        
        CALL GPXYZ(CPX,CPY,CPZ,CET,RAD,XV)
c        S_ET_TOT = S_ET_TOT + CET
      ENDDO

C ---------------
C --- ETAMAX ---
C ---------------
      ETAMAX      = NO_USE
      do i = 1, coutab(CConSa)
         call FetTab(CConSa,ID,I)
         if(abs(CConSa_z-vtx(3)).gt.0.01)then
            angle = atan2(sqrt(CConSa_x*CConSa_x+CConSa_y*CConSa_y),
     +           CConSa_z-vtx(3))
            if(angle.lt.0.) angle = angle + pi
         else
            angle = pi/2.
         endif
         eta_con = -log(tan(angle/2.))
         if(CConSa_e.gt.thres)then
            if(eta_con.gt.etamax)
     +       etamax = eta_con
         endif
      enddo 

C       Subroutine CALORTIME(Ierr)
C
C Author: Gareth Howell  - 29/5/97
C
C OUTPUT:  Ierr = 0  all fine
C          Ierr = 1  CTIME2 missing, set lateness = z vertex = 0
C **************************************************************
C +CDE, PARTAP,ZREVT.
C +CDE, CTIME2.
C +CDE, CALTuIMEC.
C        REAL Ecut,ImBcut,Vtx(3),Lateness,Avgt(0:5), AvgTErr(0:5),
C      +      TChiSq(0:5),TProb(0:5),ESum(0:5)
C        INTEGER NPMT(0:5)
C        INTEGER Ierr
C        DATA imbcut/0.35/
C        SAVE imbcut
 
 
          Ierr   = 0
          vtx(1) = 0.
          vtx(2) = 0.
c
       IF (COUTAB(CTIME2).gt.0) THEN
          CALL FETTAB(CTIME2,ID,1)
          FCAL_TM = CTime2_avtime(1) ! fcal time
          BCAL_TM = CTime2_avtime(2) ! bcal time
          RCAL_TM = CTime2_avtime(3) ! rcal time
          GLOB_TM = CTime2_avtime(4) ! global cal time
       ELSE
          FCAL_TM = -999.
          BCAL_TM = -999.
          RCAL_TM = -999.
          GLOB_TM = -999.
       ENDIF
C
      RETURN
      END
 
C     ========================================
      SUBROUTINE GPXYZ(CPX,CPY,CPZ,CET,RAD,XV)
C     ========================================
C     Auxiliary routine for angle calculation.
C
C     A.Quadt, May 1997.
C     -----------------------------------------
 
      IMPLICIT NONE
 
#include "caltru.inc"
#include "partap.inc"
 
C
C -------> Caltru is Fettab'd
C
      INTEGER    IERR
      REAL       ZP,YP,CPX,CPY,CPZ,CET,RAD,XV(3),THETA,PHI,
     1           SINPHI,SINTH,COSTH,COSPHI,XP
 
      CALL CCCXYZ (CALTRU_CELLNR,XP,YP,ZP,IERR)
      CALL CCCAPO (XP-XV(1),YP-XV(2),ZP-XV(3),RAD,THETA,PHI)
      RAD    = SQRT((XP-XV(1))**2+(YP-XV(2))**2+(ZP-XV(3))**2)
      SINPHI = SIN(PHI)
      COSPHI = COS(PHI)
      SINTH  = SIN(THETA)
      COSTH  = COS(THETA)
 
      CPX    = CALTRU_E*SINTH*COSPHI
      CPY    = CALTRU_E*SINTH*SINPHI
      CPZ    = CALTRU_E*COSTH
      CET    = CALTRU_E*SINTH
 
      RETURN
      END
 
