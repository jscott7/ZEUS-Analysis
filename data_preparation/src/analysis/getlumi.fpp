C     =================
      SUBROUTINE GETLUMI
C     =================
C     This routine reconstructs the energy and position of deposits
C     in the Lumi calorimeter (electron and photon), the BPC,
C     the 8 meter tagger and the 44 meter tagger.
C
      IMPLICIT NONE

#include "partap.inc"
#include "lmresu.inc"
#include "lmposd.inc"
#include "lmeb.inc"
c#include "t8_data.inc"
#include "common.inc"
C
      INTEGER  NO_USE, I,J,K, Ierr, Irun, JBYT
      EXTERNAL JBYT
      REAL     XPOS, ERR_X, YPOS, ERR_Y
      DATA     NO_USE / -999.99 /

C =========================
C === LUMI calorimeters ===
C =========================
      TEMPLUMG = -500
      TEMPLUME = -500
      ENE44M = 1000
      POSX_LG = NO_USE
      POSY_LG = NO_USE
      POSX_LE = NO_USE
      POSY_LE = NO_USE
 
      ENE_LG  = NO_USE
      ENE_LE  = NO_USE


C =======================
C === 44 meter tagger ===
C =======================
      TAG44_E  = 0
c      TAG44_X1 = 0
c      TAG44_X2 = 0
c      TAG44_X3 = 0
      CALL FETTAB(LMEB, ID, 1)
      TAG44_E =255-jbyt(LMEB_Rexy(8),25,8)

C -------------------
C --- LUMI ENERGY ---
C -------------------
      IF (COUTAB(LMRESU).GT.0) THEN
          CALL FETTAB(LMRESU,ID,1)
          ENE_LG = LMRESU_Enrg
          ENE_LE = LMRESU_Enre
c       else
c          ene_lg = -999
c          ene_le = -999
      ENDIF

C ----------------
C --- POSITION ---
C ----------------
	xpos = -1000.
	ypos = -1000.
      CALL LMPOSE(xpos,err_x,ypos,err_y,Ierr)
      IF (Ierr .GT. 0) THEN
         POSX_LE = XPOS
         POSY_LE = YPOS
      ENDIF
    
	xpos = -1000.
	ypos = -1000.
      CALL LMPOSG(xpos,err_x,ypos,err_y,Ierr)
      IF (Ierr .GT. 0) THEN
         POSX_LG = XPOS
         POSY_LG = YPOS
      ENDIF

C ======================
C === 8 meter tagger ===
C ======================
c      TAG8_X = NO_USE
c      TAG8_E = NO_USE
c      TAG8_P = NO_USE
c      TAG8_V = NO_USE

c      Irun = RUN_NUM
c      CALL GT_TAG8(irun,ierr)  !!! routine by Ulrike Wollmer

c      IF (Ierr.GE.0) THEN

c        TAG8_X = t8_x
c        TAG8_E = t8_e
c        TAG8_P = t8_pres
c        TAG8_V = t8_veto

c      ENDIF

C =======================
C === 44 meter tagger ===
C =======================
c      TAG44_E  = 0
c      TAG44_X1 = 0
c      TAG44_X2 = 0
c      TAG44_X3 = 0

c      IF (COUTAB(LMEB).GT.0) THEN

c      CALL FETTAB(LMEB, ID, 1)
c      TAG44_E =255-jbyt(LMEB_Rexy(8),25,8) ! ene deposit in cal

c      TAG44_X1=255-jbyt(LMEB_Rexy(8),17,8) ! ene deposit in scinti.
C                                            finger close to the bpipe

c      TAG44_X2=255-jbyt(LMEB_Rexy(8),9,8)  ! ene deposit in scinti.
C                                            finger ... next one ..

c      TAG44_X3=255-jbyt(LMEB_Rexy(8),1,8)  ! ene deposit in scinti
C                                            finger far away from bpipe
c      ENDIF
      
c      call fettab(LMEB,ID,1)
      RETURN
      END








