C     ============================
      SUBROUTINE BGDFIL(REAL_DATA)
C     ============================

      IMPLICIT NONE

#include "partap.inc"
#include "vehits.inc"
#include "common.inc"
c#include "gentup.inc"
c#include "bgdtup.inc"
C
      INTEGER  IFL,IFLX,IEL,IW,IEE
      INTEGER  Ierr
      REAL     VTX(3),ESUM
      INTEGER  HFLAG,MFLAG,AFLAG,GFLAG
      LOGICAL  REAL_DATA
      INTEGER  I, BGAS_CUT, C5(4)

C ----------------------------------------
C --- general run range classification ---
C ----------------------------------------
      CALL EVTAKE    (EVTwant)
      CALL EXOTAKE   (RUN_NUM,EVTwant,EXOwant)
      CALL XSECTAKE96(RUN_NUM,EV_NO,  XSECwant)


C --------------------
C --- QED-Comptons ---
C --------------------
      QEDC_NO = 0
      CALL SETVTX  (VTX)
      CALL COMCOS  (IFL)           !!! Tatsu T. routine.
      CALL QEDCOMP (VTX,IFLX,ESUM) !!! Thomas D. routine.

      IF (IFL.GE.1)  QEDC_NO = QEDC_NO +   1
      IF (IFL.EQ.2)  QEDC_NO = QEDC_NO +  10
      IF (IFLX.EQ.1) QEDC_NO = QEDC_NO + 100

      CALL QEDDECIDE(Ierr,VTX)
      QEDC_NO = QEDC_NO + 1000*Ierr

C --------------------------
C --- ISITAMU and MUTRIG ---
C --------------------------
      MUON_NO = 0
      MFLAG   = 0
      HFLAG   = 0
      AFLAG   = 0

      CALL ISITAMU (HFLAG)
C-------------------------------------------------------------------
C              Hflag = -1  for a non-muon                          -
C                    =  0  for unclear objects                     -
C                    =  1  something on the verge of being a muon  -
C                    =  2  potential muon                          -
C-------------------------------------------------------------------

      IF (.not.REAL_DATA) THEN
        IFL = 10
      ELSE
        IFL = 20
      ENDIF
      CALL MUTRIG(1000.0,1000.0,0.0,40.0,50.0,IFL,MFLAG)
C-------------------------------------------------------------------
C         muflag : Decision flag from mutrig, 50 for cosmic        -
C                  mu, 55 for mu through vertex, and 60 for        -
C                  beam associated halo muons.                     -
C-------------------------------------------------------------------

      CALL ALHALO2 (AFLAG)
C-------------------------------------------------------------------
C         aflag = 10 * #idendified 'halo muon modules'             -
C-------------------------------------------------------------------

      Gflag = 0
      CALL GAZZAMU(VTX,Gflag)
C---------------------------------------------------
C     GH muon finder: Gflag = 0 -> no muon found; -
C                     Gflag = 1 ->    muon found  -
C---------------------------------------------------
      MUON_NO =       1 * (HFLAG+1)
     &          +   100 *  MFLAG
     &          +  1000 *  AFLAG
     &          + 10000 *  GFLAG


C -----------------
C --- VETO WALL ---
C -----------------
      VWALL_NO = MIN(999,COUTAB(VEHITS))


C ------------------
C --- C5 COUNTER ---
C ------------------
      BGAS_CUT = 7

      CALL O1C5TM (C5(1),C5(2),C5(3),C5(4))
      C5_NO = 0
      DO I = 1,4
         IF (C5(I).GE.BGAS_CUT) C5_NO = 1
      ENDDO

      RETURN
      END


