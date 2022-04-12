C     ===================================
      SUBROUTINE GETSINISTRA(XV,YEAR_TMP)
C     ===================================
C
C -------------------------------------------------------------
C --- Run SINISTRA.                                         ---
C --- Select "best" candidate.                              ---
C --- Fill electron candidate.                              ---
C -------------------------------------------------------------

      IMPLICIT NONE

#include "partap.inc"
#include "caltru.inc"
#include "sidat95.inc"
c#include "ctrlcm.inc"
#include "vctpar.inc"
#include "vctrhl.inc"
#include "const.inc"
#include "ZufosNT.inc"
#include "common.inc"

      LOGICAL      FIRST, trk_use
      DATA         FIRST / .TRUE. /
      INTEGER      i,J,K,l,m
      INTEGER      comp_Ecand
      REAL         XV(3), CALPOS(3), SPOS(3),SRTDC(20), ECALlocal
      REAL         TRACK_P, R_DIST, trk_cone  
      REAL         cal_energy
      REAL         ESRTD, ECORR,RDIST
      REAL         EInDummy(100), ENotInDummy(100)
      INTEGER      YEAR_TMP, Ncones
      INTEGER      IDCLU
      INTEGER      IERR
      REAL         Radius(1)
      real         SINISTRACUT

      do i = 1, 56
         CellList(i) = 0
      end do

      do i = 1, 50
         do j = 1, 3
            elcells(i,j) = 0
         end do
      end do
   
      eprob = 0.
      cal_ee = 0.
      cal_xp = 0.
      cal_yp = 0.
      cal_zp = 0.
      ecelnr = 0
      trk_i = 0.
      trk_x = 0.
      trk_y = 0.
      trk_z = 0.
      trk_p = 0.
      trk_n = 0.
      trk_v = 0.
      trk_d = 0.
      trk_c = 0.
      trk_t = 0.
      srtd_x = 0.
      srtd_y = 0.
      srtd_z = 0.
      srtd_e = 0.
      srtd_ecor = 0.
      hes_x = 0.
      hes_y = 0.
      hes_e = 0.
      hes_r = 0.
      hes_f = 0.
      candid = 0

      SINISTRACUT = 0.5
      
      IF (FIRST) THEN
          FIRST = .FALSE.
          WRITE(*,*) ' *********** SINISTRA 95 ******** '
          WRITE(*,*) ' CUT VALUE ---------- ',SINISTRACUT
          WRITE(*,*) ' ******************************** '
          DO I=1,20
            SRTDC(I) = 0.0
          ENDDO
      ENDIF

      CALL SIRA95(XV,SINISTRACUT, IERR)

      IF (Ierr .NE. 0)  then
         cal_ee = -999.
         RETURN
      end if

c ----------------------------------------------
c Use Findis option 1. At low Q**2 most 
c events are by beam pipe and do not have
c tracks.
c ----------------------------------------------
      CALL FINDIS95(1, SINISTRACUT, Icand, Ierr)

      if (Ierr.ne.0) then
         cal_ee = -999.
         return
      end if
     
      comp_Ecand = Icand
      CANDID = ncand

c ---------------------------------------
c If electron found go on and get useful 
c stuff from SINISTRA, hadronic stuff from
c ZufosNT and CorandCut and SRTD info.
c ---------------------------------------
      call setvtx(XV)
      IF (ICand .GT. 0) THEN

C ---   Do tower to cell island splitting if  ---
C ---   within 60cm of RCAL beam pipe         ---
c         IF ((CANDAT(5,Icand) .LT. -140.0) .AND.
c     &        (CANDAT(3,Icand)**2 + CANDAT(4,Icand)**2) .LT.3600) THEN
            CALL Tow_to_cIsland(XV,Icand, Ierr)
c         ENDIF

C ---------------------------------
C ---   Calorimeter information ---
C ---------------------------------
         EPROB  = CANDAT(1,Icand)
         CAL_EE = CANDAT(2,Icand)
         CAL_XP = CANDAT(3,Icand)
         CAL_YP = CANDAT(4,Icand)
         CAL_ZP = CANDAT(5,Icand)
         ECELNR = POSDAT(3,Icand)

C ---   temp. storage of cell numbers for hadronic variables ---
         DO I=4,ECELNR+3
            ELCELLS(I-3,1) = POSDAT(I,Icand)
            CellList(I-3)   = POSDAT(I,Icand) !!! List for E-isol. code
         ENDDO


C ---------------------------
C ---   Do Track matching ---
C ---------------------------
         CALL TRK_ISO(INT(CANDAT(17,Icand)), trk_cone)
         TRK_I  = trk_cone
c ----------------------------
c Track x, y, z at CAL.
c ----------------------------
         TRK_X  = CANDAT(12,Icand)
         TRK_Y  = CANDAT(13,Icand)
         TRK_Z  = CANDAT(14,Icand)
c ----------------------------
c Energy of track
c ----------------------------
         TRK_P  = CANDAT(11,Icand)
         TRK_N  = CANDAT(15,Icand)
         TRK_V  = CANDAT(16,Icand)
c ----------------------------
c Distance Closest Approach
c 3D
         TRK_D  = CANDAT(18,Icand)
c 2D
         TRK_C  = CANDAT(19,Icand)
c ----------------------------
c Theta angle in radians
c ----------------------------
         TRK_T  = CANDAT(20,Icand)
c ----------------------------
c TRK X, Y, Z, P, D used later 
c for best x/y/z
C -----------------------------------
C ---   electron cone isolation ? ---
C -----------------------------------
         NrEcells   = ECELNR
         CALPOS(1)  = CAL_XP
         CALPOS(2)  = CAL_YP
         CALPOS(3)  = CAL_ZP      
        Radius(1)  = 0.8        !!! cone around electron
        Ncones     = 1          !!! test only one cone
 
        CALL IsoCones(CALPOS,   Radius, Ncones, XV,
     &                CellList, NrEcells,
     &                EInDummy, ENotInDummy,RElec,
     &                EMCenrgy,HACenrgy)
 
        EInList    = EInDummy   (1)
        ENotInList = ENotInDummy(1)

         cal_energy = cal_ee

C --------------------------
C ---   SRTD information ---
C --------------------------
         ECALlocal       = cal_ee
         CALL SRTDELEC(XV,CALPOS,ECALlocal,SRTDC,
     1        IDCLU,SPOS,ESRTD,ECORR,RDIST,Ierr)
         SRTD_X = SPOS(1)
         SRTD_Y = SPOS(2)
         SRTD_Z = SPOS(3)
         SRTD_E = ESRTD
         SRTD_ECOR = ECORR
         POSECO = FLOAT(Ierr)
     
C -------------------------
C ---   HES information ---
c X, Y from HESCLU
c -------------------------
         HES_X  = CANDAT( 8,Icand)
         HES_Y  = CANDAT( 9,Icand)  
c -------------------------
c Ehes from HESCLU
c -------------------------
         HES_E  = CANDAT( 7,Icand)
         HES_R  = CANDAT(10,Icand)
         HES_F  = CANDAT( 6,Icand)
 
      else 

         cal_ee = -999.

      end if
      
      RETURN
      END





