C     ======================================
      SUBROUTINE QEDDECIDE(QED_FLAG, VTX_IN)
C     ======================================
C     Routine that classifies events as QED compton candidates.
C     The subroutine QEDCREJ needs is called first and needs
C     to be part of the source coe used.
C
C     INPUT:   REAL    VTX_IN(3): (x,y,z) event vertex
C
C     OUTPUT:  INTEGER QED_FLAG: -1 -> strange event, check in detail !!!
C                                 0 -> not QEDC
C                                 1 -> candidate for   elastic QEDC
C                                 2 -> candidate for inelastic QEDC
C
C
C A.Quadt, 13.06.1998, version 1.1: first released version
C          29.06.1998, version 1.2: minor bug fix in cone isolation
C -------------------------------------------------------------

      IMPLICIT NONE

#include "qedtup.inc"


      INTEGER QED_FLAG, INDEX
      REAL    pi
      REAL    VTX_X, VTX_Y, VTX_Z, pos1(3),pos2(3),
     &        rad, theta,phi_ang, phi1, phi2, pt1, pt2,xx,yy,zz,
     &        vtx_in(3)

C -------------------------------------------------------------
C --- call auxiliary routine to evaluate several quantities ---
C -------------------------------------------------------------
      CALL QEDCREJ(vtx_in)

      QED_FLAG = 0
      pi = acos(-1.0)

      vtx_x =  vtx_in(1)
      vtx_y =  vtx_in(2)
      IF (vtx_in(3).gt.-999) then
        vtx_z = vtx_in(3)
      ELSE
        vtx_z = 0.0
      ENDIF

C -------------------------------------
C --- event with > 2 good electrons ---
C -------------------------------------
      IF (QNR_GOOD.GT.2) THEN
        QED_FLAG = -1
        RETURN
      ENDIF

C --------------------------
C --- rough preselection ---
C --------------------------
C      IF (empz_tot.gt.70.OR.empz_tot.lt.35) RETURN
C      IF (QCAL_EE(1).LT.9)                  RETURN
C      IF (QENotInList(1).GT.10)             RETURN

C --------------------------------------------
C --- elastic and inelastic QEDC selection ---
C --------------------------------------------
      IF (QEDCNCAND.LT.2) RETURN



C ---------------------
C --- elastic QEDCs ---
C ---------------------

C =======================
C === ge 2 candidates ===
C =======================
      IF (QEDCNCAND.LT.2) RETURN
      IF (QEProb(1).LT.0.90.OR.
     &    QEProb(2).LT.0.90) RETURN
      IF (QENotInList(1).GT.10.OR.
     &    QENotInList(2).GT.10) RETURN
      IF (QCAL_ZP(1).GT.200.AND.
     &    sqrt(QCAL_XP(1)**2+QCAL_YP(1)**2).LT.40) RETURN
      IF (QCAL_ZP(2).GT.200.AND.
     &    sqrt(QCAL_XP(2)**2+QCAL_YP(2)**2).LT.40) RETURN

C ================
C === box cuts ===
C ================
      if (qsrtd_e(1).GT.0) THEN
        pos1(1) = qsrtd_x(1)
        pos1(2) = qsrtd_y(1)
        pos1(3) = qsrtd_z(1)
      ELSE
        pos1(1) = qcal_xp(1)
        pos1(2) = qcal_yp(1)
        pos1(3) = qcal_zp(1)
      ENDIF

      if (qsrtd_e(2).GT.0) THEN
        pos2(1) = qsrtd_x(2)
        pos2(2) = qsrtd_y(2)
        pos2(3) = qsrtd_z(2)
      ELSE
        pos2(1) = qcal_xp(2)
        pos2(2) = qcal_yp(2)
        pos2(3) = qcal_zp(2)
      ENDIF

C -----------------------
C --- phi calculation ---
C -----------------------
      XX = pos1(1) - VTX_X !!! x-vertex
      YY = pos1(2) - VTX_Y !!! y-vertex
      ZZ = pos1(3) - VTX_Z !!! z-vertex
      CALL CCCAPO(XX,YY,ZZ,RAD,THETA,PHI_ANG)
      phi1 = phi_ang
      PT1  = QCAL_EE(1) * SIN(THETA)
      IF (phi1.GT.2.0*pi) phi1 = phi1 -2.0*pi
      IF (phi1.LT.0.0)    phi1 = phi1 +2.0*pi

      XX = pos2(1) - VTX_X !!! x-vertex
      YY = pos2(2) - VTX_Y !!! y-vertex
      ZZ = pos2(3) - VTX_Z !!! z-vertex
      CALL CCCAPO(XX,YY,ZZ,RAD,THETA,PHI_ANG)
      phi2 = phi_ang
      PT2  = QCAL_EE(2) * SIN(THETA)
      IF (phi2.GT.2.0*pi) phi2 = phi2 -2.0*pi
      IF (phi2.LT.0.0)    phi2 = phi2 +2.0*pi



C ==============================
C === EMC/ETOT fraction cuts ===
C ==============================
      IF (EEMC_TOT+EHAD_TOT           .LE.0.0) RETURN

      IF (EEMC_TOT/(EEMC_TOT+EHAD_TOT).LT.0.9) THEN
        GOTO 500
      ENDIF


C ===========================
C === ETOT-Ee-Egamma cuts ===
C ===========================
      IF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.200.OR.
     &    QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.200) THEN
        IF (eemc_tot+ehad_tot-qcal_ee(1)-qcal_ee(2).GT.3.0) GOTO 500
      ELSEIF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.-140.OR.
     &        QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.-140) THEN
        IF (eemc_tot+ehad_tot-qcal_ee(1)-qcal_ee(2).GT.3.0) GOTO 500
      ELSE
        IF (eemc_tot+ehad_tot-qcal_ee(1)-qcal_ee(2).GT.3.0) GOTO 500
      ENDIF

C =================
C === dphi cuts ===
C =================
      IF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.200.OR.
     &    QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.200) THEN
        IF (abs(phi1-phi2).GT.3.50.OR.
     &      abs(phi1-phi2).LT.2.78) GOTO 500
      ELSEIF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.-140.OR.
     &        QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.-140) THEN
        IF (abs(phi1-phi2).GT.3.30.OR.
     &      abs(phi1-phi2).LT.2.90) GOTO 500
      ELSE
        IF (abs(phi1-phi2).GT.3.80.OR.
     &      abs(phi1-phi2).LT.2.40) GOTO 500
        IF ((abs(phi1-phi2).GT.3.30.OR.
     &       abs(phi1-phi2).LT.2.90).AND.
     &       QCAL_ZP(1).LT.-140.AND.
     &       QCAL_ZP(2).LT.-140.AND.
     &       sqrt(QCAL_xp(1)**2+QCAL_YP(1)**2).LT.50.AND.
     &       sqrt(QCAL_XP(2)**2+QCAL_YP(2)**2).LT.50) THEN
          IF (ABS(pt1-pt2).GT.3.0) GOTO 500
        ENDIF
      ENDIF

      IF (QCAL_EE(2).LT.7.OR.
     &    QCAL_ee(1).LT.10) THEN
        IF (ABS(pt1-pt2).GT.3.0) GOTO 500
      ENDIF

C =============================
C === set elastic QEDC flag ===
C =============================
      QED_FLAG = 1
      RETURN




 500  CONTINUE
C -----------------------------
C --- check inelastic QEDCs ---
C -----------------------------

C --------------------
C --- delta E cuts ---
C --------------------

C --- FCAL object ---
      IF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.200.OR.
     &    QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.200) THEN
          IF (QCAL_ZP(1).GT.200) THEN
            index = 1
          ELSE
            index = 2
          ENDIF
          rad = sqrt(QCAL_xp(index)**2+qcal_yp(index)**2)
          IF (rad.GT.90) THEN
            IF (eemc_tot+ehad_tot
     &          -qcal_ee(1)-qcal_ee(2)-FCALE90.GT.3.0) RETURN
          ELSEIF (rad.GT.70) THEN
            IF (eemc_tot+ehad_tot
     &          -qcal_ee(1)-qcal_ee(2)-FCALE70.GT.3.0) RETURN
          ELSEIF (rad.GT.50) THEN
            IF (eemc_tot+ehad_tot
     &          -qcal_ee(1)-qcal_ee(2)-FCALE50.GT.3.0) RETURN
          ELSE
C ---       currently no cut for FCAL objext within 50 cm
            RETURN
          ENDIF

C --- BCAL object ---
      ELSEIF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.-140.OR.
     &        QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.-140) THEN
        IF (eemc_tot+ehad_tot
     &      -qcal_ee(1)-qcal_ee(2)-FCALEm.GT.3.0) RETURN
C --- RCAL object ---
      ELSE
        IF (eemc_tot+ehad_tot
     &      -qcal_ee(1)-qcal_ee(2)-FCALEm.GT.5.0) RETURN
      ENDIF


C ===========================
C === energy fraction cut ===
C ===========================
C --- FCAL objext ---
      IF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.200.OR.
     &    QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.200) THEN
          IF (QCAL_ZP(1).GT.200) THEN
            index = 1
          ELSE
            index = 2
          ENDIF
          rad = sqrt(QCAL_xp(index)**2+qcal_yp(index)**2)
          IF (rad.GT.90) THEN
            IF (EEMC_TOT+EHAD_TOT-fcalem.GT.0.0) THEN
             IF ((eemc_tot-fcalee90)/ 
     &           (eemc_tot+ehad_tot-fcale90).LT.0.9) RETURN
            ENDIF
          ELSEIF (rad.GT.70) THEN
            IF (EEMC_TOT+EHAD_TOT-fcalem.GT.0.0) THEN
             IF ((eemc_tot-fcalee70)/ 
     &           (eemc_tot+ehad_tot-fcale70).LT.0.9) RETURN
            ENDIF
          ELSEIF (rad.GT.50) THEN
            IF (EEMC_TOT+EHAD_TOT-fcalem.GT.0.0) THEN
             IF ((eemc_tot-fcalee50)/ 
     &           (eemc_tot+ehad_tot-fcale50).LT.0.9) RETURN
            ENDIF
          ELSE
            RETURN
          ENDIF

C --- BCAL object ---
      ELSEIF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.-140.OR.
     &        QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.-140) THEN
          IF (EEMC_TOT+EHAD_TOT-fcalem.GT.0.0) THEN
           IF ((eemc_tot-fcaleem)/
     &         (eemc_tot+ehad_tot-fcalem).LT.0.9) RETURN
          ENDIF
C --- RCAL object ---
      ELSE
          IF (EEMC_TOT+EHAD_TOT-fcalem.GT.0.0) THEN
           IF ((eemc_tot-fcaleem)/
     &         (eemc_tot+ehad_tot-fcalem).LT.0.9) RETURN
          ENDIF
      ENDIF


C ===============
C === dpt cut ===
C ===============
C --- FCAL object ---
      IF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.200.OR.
     &    QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.200) THEN
          IF (abs(pt1-pt2).GT.5) RETURN
C --- BCAL object ---
      ELSEIF (QCAL_EE(1).GT.0.0.AND.QCAL_ZP(1).GT.-140.OR.
     &        QCAL_EE(2).GT.0.0.AND.QCAL_ZP(2).GT.-140) THEN
          IF (abs(pt1-pt2).GT.7) RETURN
          IF (QCAL_EE(2).LT.7.OR.
     &        QCAL_EE(1).LT.10) THEN
            IF (ABS(pt1-pt2).GT.4.0) RETURN
          ENDIF
C --- RCAL object ---
      ELSE
          IF (abs(pt1-pt2).GT.5) RETURN
          IF (QCAL_EE(2).LT.7.OR.
     &        QCAL_EE(1).LT.10) THEN
            IF (ABS(pt1-pt2).GT.4.0) RETURN
          ENDIF
      ENDIF

C ===============================
C === check outer FCAL energy ===
C ===============================
      IF (FCALEM-FCALE50.GT.5.0) RETURN

C ===============================
C === set inelastic QEDC flag ===
C ===============================
      QED_FLAG = 2

      RETURN
      END





C     =======================
      SUBROUTINE QEDCREJ(VTX)
C     =======================
C     Auxiliary routine for QED Compton event rejection.
C
C
C A.Quadt, 10.06.1998, version 1.1: first released version
C          29.06.1998, version 1.2: minor bug fix in cone isolation
C -------------------------------------------------------------

      IMPLICIT NONE

#include "partap.inc"
#include "sidat95.inc"
#include "caltru.inc"
#include "qedtup.inc"


      INTEGER    IERR, I, J, vector(5), Icand, Iloop,
     &           Icand1, CellList(100), nr_cells, Ncones, IDCLU,
     &           NrEcells
      REAL       SINISCUT, VTX(3), PX_TOT, PY_TOT, XP, YP, ZP,
     &           max_prob, trk_cone,Radius(1), 
     &           EInDummy(100), ENotInDummy(100), RElec,dummy1, dummy2,
     &           ECAL,SPOS(3),SRTDC(20),ESRTD, ECORR, RDIST,
     &           cpx,cpy,cpz,cet,rad,calpos(3)
      PARAMETER (SINISCUT = 0.8)
      LOGICAL    FIRST
      DATA       FIRST/.TRUE./
      CHARACTER KIND*5



      IF (FIRST) THEN
         FIRST = .FALSE.
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' +-------------------------+'
         WRITE(*,*) ' | This is QEDCREJ v1.0    |'
         WRITE(*,*) ' |                         |'
         WRITE(*,*) ' | questions & comments to |'
         WRITE(*,*) ' | quadt@desy.de           |'
         WRITE(*,*) ' +-------------------------+'
         WRITE(*,*) ' '
         DO I=1,20
            SRTDC(I) = 0.0
         ENDDO
         WRITE(*,*) ' '
         WRITE(*,*) ' '
      ENDIF


      EMPZ_TOT  = 0.0
      EEMC_TOT  = 0.0
      EHAD_TOT  = 0.0
      PT_TOT    = 0.0
      PX_TOT    = 0.0
      PY_TOT    = 0.0
      FCALE50   = 0.0
      FCALE70   = 0.0
      FCALE90   = 0.0
      FCALEm    = 0.0
      FCALEE50  = 0.0
      FCALEE70  = 0.0
      FCALEE90  = 0.0
      FCALEEm   = 0.0

      DO I=1,COUTAB(CALTRU)
        CALL FETTAB(CALTRU, ID, I)
        CALL CCWHAT(CALTRU_CELLNR, KIND, VECTOR, IERR)
        CALL MY_GPXYZ(CPX,CPY,CPZ,CET,RAD,VTX)      
        EMPZ_TOT  = EMPZ_TOT + CALTRU_E - CPZ
        PX_TOT    = PX_TOT   + CPX
        PY_TOT    = PY_TOT   + CPY
        IF (KIND(1:2).EQ.'RE'.OR.
     &      KIND(1:2).EQ.'BE'.OR.
     &      KIND(1:2).EQ.'FE'.OR.
     &      KIND(1:5).EQ.'RHAC0'.OR.
     &      KIND(1:5).EQ.'FHAC0') THEN
          EEMC_TOT = EEMC_TOT + CALTRU_E
        ELSE
          EHAD_TOT = EHAD_TOT + CALTRU_E
        ENDIF
        CALL CCCXYZ (CALTRU_CELLNR,XP,YP,ZP,IERR)
        IF (KIND(1:1).EQ.'F') THEN
          IF (sqrt(xp**2+yp**2).LE.50) THEN
            FCALE50 = FCALE50 + CALTRU_E
          ENDIF
          IF (sqrt(xp**2+yp**2).LE.70) THEN
            FCALE70 = FCALE70 + CALTRU_E
          ENDIF
          IF (sqrt(xp**2+yp**2).LE.90) THEN
            FCALE90 = FCALE90 + CALTRU_E
          ENDIF
          FCALEm = FCALEm + CALTRU_E
        ENDIF
        IF (KIND(1:2).EQ.'FE'.OR.
     &      KIND(1:5).EQ.'FHAC0') THEN
          IF (sqrt(xp**2+yp**2).LE.50) THEN
            FCALEE50 = FCALEE50 + CALTRU_E
          ENDIF
          IF (sqrt(xp**2+yp**2).LE.70) THEN
            FCALEE70 = FCALEE70 + CALTRU_E
          ENDIF
          IF (sqrt(xp**2+yp**2).LE.90) THEN
            FCALEE90 = FCALEE90 + CALTRU_E
          ENDIF
          FCALEEm = FCALEEm + CALTRU_E
        ENDIF
      ENDDO

      PT_TOT = sqrt(PX_TOT**2+PY_TOT**2)


      QEDCNCAND = 0
      CALL SIRA95(VTX,SINISCUT,ierr)
      IF (Ierr .NE. 0)  RETURN
      QEDCNCAND = Ncand

      DO Iloop = 1,2
        QCAL_EE(Iloop) = 0.0
        QCAL_XP(Iloop) = 0.0
        QCAL_YP(Iloop) = 0.0
        QCAL_ZP(Iloop) = 0.0
        QSRTD_E(Iloop) = 0.0
        QSRTD_X(Iloop) = 0.0
        QSRTD_Y(Iloop) = 0.0
        QSRTD_Z(Iloop) = 0.0
        IF (Iloop.eq.1) THEN
          CALL FINDIS96( 5, SINISCUT, Icand, Ierr, VTX(3))
          Icand1 = Icand
        ELSE
          Icand    = 0
          max_prob = -999.99
          DO J=1,Ncand
            IF (J.NE.Icand1.AND.CANDAT(1,J).GT.max_prob) THEN
              max_prob = CANDAT(1,J)
              Icand    = J
            ENDIF
          ENDDO
        ENDIF


        IF (ICAND.GT.0) THEN


C ---   Do tower to cell island splitting if  ---
C ---   within 60cm of RCAL beam pipe         ---
        IF ((CANDAT(5,Icand) .LT. -140.0) .AND.
     &      (CANDAT(3,Icand)**2 + CANDAT(4,Icand)**2) .LT.3600) THEN
          CALL Tow_to_cIsland(VTX,Icand, Ierr)
        ENDIF

C ---------------------------------
C ---   Calorimeter information ---
C ---------------------------------
        QEPROB (Iloop) = CANDAT(1,Icand)
        QCAL_EE(Iloop) = CANDAT(2,Icand)
        QCAL_XP(Iloop) = CANDAT(3,Icand)
        QCAL_YP(Iloop) = CANDAT(4,Icand)
        QCAL_ZP(Iloop) = CANDAT(5,Icand)
        QECELNR(Iloop) = MIN(POSDAT(3,Icand),50)


C ---   temp. storage of cell numbers for hadronic variables ---
        DO I=1,QECELNR(Iloop)
          CellList(I) = POSDAT(3+I,Icand)   !!! List for E-isol. code
        ENDDO

       
C ---------------------------
C ---   Do Track matching ---
C ---------------------------
C        CALL TRK_ISO(INT(CANDAT(17,Icand)), trk_cone)
C        QTRK_I(Iloop)  = trk_cone
        QTRK_I(Iloop)  = 0.0
        QTRK_X(Iloop)  = CANDAT(12,Icand)
        QTRK_Y(Iloop)  = CANDAT(13,Icand)
        QTRK_Z(Iloop)  = CANDAT(14,Icand)
        QTRK_P(Iloop)  = CANDAT(11,Icand)
        QTRK_N(Iloop)  = CANDAT(15,Icand)
        QTRK_V(Iloop)  = CANDAT(16,Icand)
        QTRK_D(Iloop)  = CANDAT(18,Icand)
        QTRK_C(Iloop)  = CANDAT(19,Icand)
        QTRK_T(Iloop)  = CANDAT(20,Icand)

C -----------------------------------
C ---   electron cone isolation ? ---
C -----------------------------------
        NrEcells   = QECELNR(Iloop)
        IF (QTRK_D (Iloop).GE.0.0    .and.
     &      QCAL_ZP(Iloop).GT.-130.0 .and.
     &      QCAL_ZP(Iloop).LT. 200.0) THEN
          CALPOS(1)  = QTRK_X(Iloop)
          CALPOS(2)  = QTRK_Y(Iloop)
          CALPOS(3)  = QTRK_Z(Iloop)
        ELSE
          CALPOS(1)  = QCAL_XP(Iloop)
          CALPOS(2)  = QCAL_YP(Iloop)
          CALPOS(3)  = QCAL_ZP(Iloop)
        ENDIF
        Radius(1)  = 0.8        !!! cone around electron
        Ncones     = 1          !!! test only one cone

        CALL IsoCones(CALPOS,   Radius, Ncones, VTX,
     &                CellList, NrEcells,
     &                EInDummy, ENotInDummy,RElec,
     &                dummy1, dummy2)

        QEInList   (Iloop) = EInDummy   (1)
        QENotInList(Iloop) = ENotInDummy(1)


        CALPOS(1)  = QCAL_XP(Iloop)
        CALPOS(2)  = QCAL_YP(Iloop)
        CALPOS(3)  = QCAL_ZP(Iloop)

C --------------------------
C ---   SRTD information ---
C --------------------------
        ECAL       = QCAL_EE(Iloop)
        CALL SRTDELEC(VTX,CALPOS,ECAL,SRTDC,
     1                IDCLU,SPOS,ESRTD,ECORR,RDIST,Ierr)
        QSRTD_X(Iloop) = SPOS(1)
        QSRTD_Y(Iloop) = SPOS(2)
        QSRTD_Z(Iloop) = SPOS(3)
        QSRTD_E(Iloop) = ESRTD
        QPOSECO(Iloop) = FLOAT(Ierr)

        ENDIF
      ENDDO

      qnr_good = 0
      DO J=1,Ncand
        IF (CANDAT(1,J).GT.0.9) THEN
          CALPOS(1)  = CANDAT(3,J)
          CALPOS(2)  = CANDAT(4,J)
          CALPOS(3)  = CANDAT(5,J)
          NrEcells   = MIN(POSDAT(3,J),50)

C ---     temp. storage of cell numbers for hadronic variables ---
          DO I=1,NrEcells
            CellList(I) = POSDAT(3+I,J)   !!! List for E-isol. code
          ENDDO
          Radius(1)  = 0.8        !!! cone around electron
          Ncones     = 1          !!! test only one cone

          CALL IsoCones(CALPOS,   Radius, Ncones, VTX,
     &                  CellList, NrEcells,
     &                  EInDummy, ENotInDummy,RElec,
     &                  dummy1, dummy2)
          IF (ENotInDummy(1).LE.5.0) THEN
            qnr_good = qnr_good + 1
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END

