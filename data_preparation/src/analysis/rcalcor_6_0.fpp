C
C Version 1 Author: N. Tuning
C          6/4/98
C
C Version 2  Modification: M. Wodarczyk 9-June-98
C       Modifiy cells with the DA method with improved statistics.
C       Changed overall RCAL scale to 2.2 as seen from outer RCAL.
C
C Version 3 Modification: M. Wodarczyk 13-Aug-98
C       Corrected inner rcal cells from Niles plots.  
C       Changes from previous cell values:
C         1. New box cut "H1 cut" to remove beampipe dead material problem.
C         2. Low Mip KP evemts used.
C
C Version 4.0 Modification: M. Wodarczyk 2-Sep-98
C       Cell 37538 was mis-entered by hand in version 3.0.
C         the wrong correction of 4.819% was corrected to 2.45%.
C       Version information moved to top of program for easy access.
C
C Version 5.0 Modification: M. Wodarczyk 24-Sep-98
C       1) Make Bcal correction 5% as in Calcorr
C       2) Make routine compatable with EM by calling
C            emSetCaltruFactors.  This means that to compile
C            you need at least emcaltruf.fpp and empar.inc
C            They are in em version 1.4 or in
C /afs/desy.de/group/zeus.zsmsm/ZEUSSysSoft/Released/zeus/ZES/zesphys/v1998a/src/fsig/em/emcaltruf.fpp
C /afs/desy.de/group/zeus.zsmsm/ZEUSSysSoft/Released/zeus/ZES/zesphys/v1998a/inc/fsig/empar.inc
C
C 
C Version 6.0 Modification: M. Wodarczyk 09-Feb-99
C       1) corrected for new crack positions as measured by Rolf Deffner in
C          the new non-uniformity correction
C       2) Corrected cells near Rcal BP where the old non-uniformity and
C          Rcalcorr double corrected.  
C       3) Enlarged Box cut near rcal BP so KP electrons leaking into the 
C          BP do not produce a large correction.  Some sort of additional 
C          leakage correction will need to be applied later based on 
C          the electron position and energy.
C
C
C


C     ===========================
      SUBROUTINE RcalCorr6(status)
C     ===========================
c     Modification of calCorr. 
c     Scale factors for individual RCAL EMC Cells are implemented 
c     ====> See routine RcalCor for documentation.
c 
c     NOTE: Calls Subroutine RcalCor( Poser, Factor)
c           input:  Poser
c           output: Factor
c     *****************************************************************
c       Author:  N. Tuning
c       Version: 22/04/98
c
c       Modification: M. Wodarczyk 9-June-98
c         Make clear that monte carlo is not corrected
c
c     *****************************************************************
c
c
C     ==========================
C      subroutine calCorr(status)
C     ==========================y
C Apply 1.06 to BCAL energy and 1.025 to RCAL energy
C
c     *****************************************************************
c       Author:  Jon Labs
c       Version: 1.0 (22/05/96)
c
c       Purpose: Apply corrections to cal energy in the data
c                to account for the 6% discrepency observed in BCAL and
c                2.5% in RCAL in the absolute energy scale
C**********************************************************************
C
C                for 95 data and for MC>=num95v2.1
C          the present recommendation is to apply 5% in the BCAL
C
C**********************************************************************
C 1) SRTDELEC, or any other existing version of SRTD energy correction, should
C   be called before CALCORR, otherwise the electron energy returned by
C   SRTDELEC will be over corrected.  Wouter Verkerke will provide a
C   modification to his routine soon so that this will no longer be an issue.
C 2) It cannot be called before corrections based on KP data, for similar
C    reasons.
C 3) Noisup94, or any existing noise suppression routine, should be called
C   before CALCORR.  A volunteer is needed to study the effect on Noisup94 when
C   CALCORR is called first.
C 4) Sinistra, and other electron finders, should be called before CALCORR.
C   The effect of CALCORR on the electron finders, at least Sinistra, needs to
C   be studied, and we are looking for a volunteer.
C 5) The island finder should be called before CALCORR.  A volunteer
C   is needed to study the effect of calling CALCORR before the island finder.
C 6) The condensates table, CCONSA, would have to be rebuilt.
C
C
C If the routine is called twice for the same event, it will not apply the
C correction the second time.  This routine will work on the RDST's and
C miniDST's, and the tables it requires are of course the CALTRU tables and the
C CCOR table.
C************************************************************************************
c
c       Input:  none
c       Output:
c                status =  1 if routine was called and executed successfully
c                            but no correction was applied
c                          0 if routine was called and executed successfully
c                            and a correction was applied
c                         -1 if there is no record as to whether or not
c                            the correction has already been made; no
c                            correction is then made
c
c
c     *****************************************************************
 
      implicit none
 
#include "partap.inc"
#include "zdskey.inc"
#include "caltru.inc"
#include "fmckin.inc"
#include "ccor.inc"
 
      integer status

c
c local variables
c
      integer cnt
      Logical MonteCarlo
      real scale

c
c give the defaults return status
c
      status = 0
 
c
c I don't want to allow the user to call both calcorr and rcalcorr 
c on the same event.  calcorr needs CCOR to exist and puts a value
c in ccor_foffset that is not 100.
c if ccor_foffset == 100 then I am safe that calcorr was not called.
c
      if (coutab(CCOR).le.0) then
         status = -1.
         return
      endif
 
      call fettab(CCOR, id, 1)
      if (ccor_foffset.ne.100.) then
         status = 1
         return
      endif

c
c now record that this correction has been made so that it is not
c made again.
c
      ccor_foffset = 99999  ! mark ccor_foffset as changed (i.e. not 100)
      call reptab(CCOR)
 
 
      MonteCarlo = (coutab(FMCKIN).gt.0)
 
c
c okay at this point we should be looking at data without the correction
c start looping through caltru.
c
      if (.not.MonteCarlo) then
         cnt = coutab(CALTRU)
         do while(cnt.gt.0)
            
            call fettab(caltru, id, cnt)
            
            
            CALL RCALCOR( CALTRU_CELLNR, SCALE)
            caltru_e = caltru_e*scale
            caltru_imbal = caltru_imbal*scale
            
            call reptab(CALTRU)
            
 10         continue
            
            cnt = cnt - 1
            
         enddo
      endif  ! not MonteCarlo
 
      return
      end
 
 
 
       
      Subroutine RCALCOR( POSER, FACTOR)
      
C  -----------------------------------------------------
C  Input:  Cell number.
C  Output: Correction factor.
C
C  This routine provides correction factors for
C  individual RCAL cells.
C  The factors are obtained in two ways:
C
C  1) Kinematic Peak events.
C     The peak position in data and MC is fitted 
C     separately. From the difference the correction
C     factor is calculated.
C     The 44 inner cells use the KP correction factors.
C     (precision within 0.5 %)
C
C  2) Double Angle method.
C     The energy of the electron is predicted by the 
C     da-formula. The difference between measured and 
C     predicted (corrected for true/predicted) gives
C     the correction factor.
C     Cells without KP information use the DA method,
C     70 cells are corrected.(precision within 1 - 1.5%)
C
C Note: Only electrons in the center of the cell were
C       taken into account. (2 cm in x and 1 cm in y)
C       The electron energy is corrected for dead
C       material (72 MeV/mip in PRES,28 MeV/mip in SRTD)
C       both for data and for Monte Carlo.
C
C  If no correction factor for the cell was obtained
C  (if low statistics in RCAL or if B/FCAL cell), the
C  following correction factors are given:
C
C  RCAL: 1.022    (The mean of the obtained factors !)
C  BCAL: 1.050    (See CALCOR)
C  FCAL: 1.000    (See CALCOR)
C
C  See top of file for version information.
C
C -----------------------------------------------------
       
        IMPLICIT NONE
       
        INTEGER I
        INTEGER POSER
        REAL    FACTOR
        REAL    CELL_INFO( 2,  124)
        INTEGER NCELLS
        DATA    NCELLS/  124 /
        logical first
        Save First, cell_info
        Data First /.True./
 
        REAL dummy1(2,25)
        DATA dummy1/
     &       36498 ,  0.9740 , 
     &       36500 ,  0.9832 , 
     &       36514 ,  0.9930 , 
     &       36516 ,  1.0012 , 
     &       36530 ,  0.9719 , 
     &       36532 ,  1.0066 , 
     &       36546 ,  1.0111 , 
     &       36562 ,  0.9753 , 
     &       36996 ,  0.9897 , 
     &       37010 ,  0.9938 , 
     &       37012 ,  1.0036 , 
     &       37026 ,  1.0248 , 
     &       37028 ,  0.9779 , 
     &       37042 ,  1.0131 , 
     &       37044 ,  0.9725 , 
     &       37058 ,  0.9837 , 
     &       37060 ,  0.9847 , 
     &       37074 ,  1.0034 , 
     &       37076 ,  0.9804 , 
     &       37090 ,  0.9633 , 
     &       37492 ,  0.9939 , 
     &       37506 ,  1.0180 , 
     &       37508 ,  1.0149 , 
     &       37522 ,  1.0260 , 
     &       37524 ,  0.9945 / 
        REAL dummy2(2,25)
        DATA dummy2/
     &       37538 ,  1.0161 , 
     &       37540 ,  1.0123 , 
     &       37554 ,  1.0178 , 
     &       37556 ,  0.9927 , 
     &       37570 ,  1.0439 , 
     &       37572 ,  0.9961 , 
     &       37586 ,  1.0007 , 
     &       37588 ,  0.9739 , 
     &       37602 ,  0.9879 , 
     &       37604 ,  1.0748 , 
     &       37618 ,  1.0175 , 
     &       37620 ,  1.0297 , 
     &       38002 ,  1.0328 , 
     &       38004 ,  1.0361 , 
     &       38018 ,  1.0632 , 
     &       38020 ,  1.0100 , 
     &       38034 ,  1.0190 , 
     &       38036 ,  0.9957 , 
     &       38050 ,  1.0291 , 
     &       38052 ,  1.0292 , 
     &       38066 ,  0.99642 , 
     &       38068 ,  1.00396 , 
     &       38082 ,  1.00954 , 
     &       38084 ,  1.0560 , 
     &       38098 ,  1.0403 / 
        REAL dummy3(2,25)
        DATA dummy3/
     &       38100 ,  1.0055 , 
     &       38114 ,  1.0099 , 
     &       38116 ,  1.0082 , 
     &       38130 ,  1.0189 , 
     &       38514 ,  1.0325 , 
     &       38516 ,  1.0306 , 
     &       38530 ,  1.0362 , 
     &       38532 ,  1.0307 , 
     &       38546 ,  1.0348 , 
     &       38548 ,  1.0226 , 
     &       38562 ,  1.0135 , 
     &       38564 ,  1.0109 , 
     &       38594 ,  0.99806 , 
     &       38596 ,  1.0459 , 
     &       38610 ,  1.0225 , 
     &       38612 ,  1.0482 , 
     &       38626 ,  1.0381 , 
     &       38628 ,  1.0577 , 
     &       38642 ,  1.0097 , 
     &       38644 ,  1.0320 , 
     &       39026 ,  1.0220 , 
     &       39028 ,  1.0406 , 
     &       39042 ,  1.0134 , 
     &       39044 ,  1.0485 , 
     &       39058 ,  1.1322 / 
        REAL dummy4(2,25)
        DATA dummy4/
     &       39060 ,  1.0209 , 
     &       39074 ,  1.0302 , 
     &       39076 ,  1.0346 , 
     &       39090 ,  1.0385 , 
     &       39092 ,  1.00722 , 
     &       39106 ,  1.01948 , 
     &       39108 ,  1.0424 , 
     &       39122 ,  1.0313 , 
     &       39124 ,  1.0193 , 
     &       39138 ,  1.0045 , 
     &       39140 ,  1.0432 , 
     &       39154 ,  1.0045 , 
     &       39540 ,  0.9797 , 
     &       39554 ,  0.9796 , 
     &       39556 ,  1.0119 , 
     &       39570 ,  1.0063 , 
     &       39572 ,  1.0099 , 
     &       39586 ,  0.9709 , 
     &       39588 ,  0.9540 , 
     &       39602 ,  0.9928 , 
     &       39604 ,  0.9859 , 
     &       39618 ,  0.9942 , 
     &       39620 ,  0.9963 , 
     &       39634 ,  1.0042 , 
     &       39636 ,  0.9831 / 
        REAL dummy5(2,24)
        DATA dummy5/
     &       39650 ,  1.0139 , 
     &       39652 ,  1.0442 , 
     &       39666 ,  0.9830 , 
     &       40066 ,  0.9839 , 
     &       40068 ,  1.0438 , 
     &       40082 ,  1.0083 , 
     &       40084 ,  1.0053 , 
     &       40098 ,  1.0357 , 
     &       40100 ,  1.0169 , 
     &       40114 ,  1.0097 , 
     &       40116 ,  1.0206 , 
     &       40130 ,  1.0238 , 
     &       40132 ,  1.0609 , 
     &       40146 ,  1.0198 , 
     &       40148 ,  0.9877 , 
     &       40162 ,  1.0277 , 
     &       40164 ,  1.0220 , 
     &       40610 ,  0.9899 , 
     &       40612 ,  0.9950 , 
     &       40626 ,  0.9539 , 
     &       40628 ,  1.0123 , 
     &       40642 ,  0.9617 , 
     &       40644 ,  0.9949 , 
     &       40660 ,  1.0390 / 
       

      if (first) then
         first = .FALSE.

C
C  Tell EM what factor I use for global corrections below
C
         call emSetCaltruFactors(1.00, 1.05, 1.022)

         DO I = 1, 25
            CELL_INFO( 1, I+ 0) = dummy1(1,I)
            CELL_INFO( 2, I+ 0) = dummy1(2,I)
            CELL_INFO( 1, I+25) = dummy2(1,I)
            CELL_INFO( 2, I+25) = dummy2(2,I)
            CELL_INFO( 1, I+50) = dummy3(1,I)
            CELL_INFO( 2, I+50) = dummy3(2,I)
            CELL_INFO( 1, I+75) = dummy4(1,I)
            CELL_INFO( 2, I+75) = dummy4(2,I)
         ENDDO
         DO I = 1,  24
            CELL_INFO( 1, I+100) = dummy5(1,I)
            CELL_INFO( 2, I+100) = dummy5(2,I)
         ENDDO
      endif
       
       IF (POSER .LT. 16384) THEN       ! FCAL
            FACTOR = 1.000
       ELSEIF (POSER .LT. 32768) THEN   ! BCAL
            FACTOR = 1.050
       ELSEIF (POSER .GE. 32768) THEN   ! RCAL
            FACTOR = 1.022
            DO I=1, NCELLS
              IF (POSER .EQ. CELL_INFO(1, I)) THEN
                 FACTOR = CELL_INFO(2, I)
              ENDIF
            ENDDO
       ELSE
            FACTOR = 1.000
       ENDIF
       
       
       RETURN
       END
       
