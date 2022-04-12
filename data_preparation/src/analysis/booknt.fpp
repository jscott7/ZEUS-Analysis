C     ============================
      SUBROUTINE BOOKNT
C     ============================
C     In this routine the NTUPLE is booked. If Running on MC
C     the block MR1TUP is replace by MR2tup. By commenting out
C     certain declaration commands the size and the entries in
C     the NTUPLE can be chosen.
C
C     A.Quadt, May 1997.
C     ---------------------------------------------------------

      IMPLICIT NONE

#include "ezhbook.inc"
#include "common.inc"

      INTEGER Ierr, LENOCC
      INTEGER IQUEST
      COMMON/QUEST/IQUEST(100)

      IQUEST(10)=65000

      CALL HLIMIT(-Nwds_HBOOK)
      CALL HROPEN(80,'ntuple','ntuple.rz','N', 4096, Ierr)
      IF (Ierr.NE.0) THEN
        WRITE(*,*) 'Error in opening NTUPLE !'
        STOP
      ENDIF

      CALL HBNT(10, 'ntuple', 'D')
      CALL HBNAME(10,flt_form(1),TrigDat, flt_form(2))
      CALL HBNAME(10,tlt_form(1),TLT, tlt_form(2))
      CALL HBNAME(10,trk_form(1),VCT_XVC, trk_form(2))
      CALL HBNAME(10,cal_form(1),FEMC_EN, cal_form(2))  
      CALL HBNAME(10,elec_form(1),EPROB, elec_form(2))
      CALL HBNAME(10,zufos1_form(1),ZPxBSp, zufos1_form(2))
      CALL HBNAME(10,zufos2_form(1),Zgamma, zufos2_form(2))
      CALL HBNAME(10,zufos3_form(1),ZEminPz, zufos3_form(2))
      CALL HBNAME(10,zufos4_form(1),ZPxCal, zufos4_form(2))
      CALL HBNAME(10,zufos5_form(1),ZufoE, zufos5_form(2))
      CALL HBNAME(10,temp_form(1),TEMPLUMG, temp_form(2))
      CALL HBNAME(10,tag1_form(1),ENE44M, tag1_form(2))
      CALL HBNAME(10,tag2_form(1),TAG44_E, tag2_form(2))
      CALL HBNAME(10,lumi_form(1),ENE_LE, lumi_form(2))
      CALL HBNAME(10,lumi2_form(1),ENE_LG, lumi2_form(2))
      CALL HBNAME(10,gen_form(1),RUN_NUM, gen_form(2))
      CALL HBNAME(10,mctrue_form(1), MCELE, mctrue_form(2))
      call hbname(10,BGDTUP_form(1),EVTwant,BGDTUP_form(2))

      call HPRNT(10)
      WRITE(*,*)
      WRITE(*,*) '  NTUPLE succesfully booked ...'
      WRITE(*,*)

      RETURN
      END





