      INTEGER    MAXQ
      PARAMETER (MAXQ = 2)
      COMMON / QEDTUP / QEDCNCAND,
     &                  EMPZ_TOT,EEMC_TOT,EHAD_TOT,PT_TOT,
     &                  FCALE50,  FCALE70,  FCALE90, FCALEm,
     &                  FCALEE50, FCALEE70, FCALEE90,FCALEEm,
     &                  QEPROB(MAXQ),     QECELNR(MAXQ),
     &                  QCAL_EE(MAXQ),    QCAL_XP(MAXQ), QCAL_YP(MAXQ),
     &                  QCAL_ZP(MAXQ),
     &                  QEInList(MAXQ),   QENotInList(MAXQ),
     &                  QSRTD_E(MAXQ),    QSRTD_X(MAXQ),
     &                  QSRTD_Y(MAXQ),
     &                  QSRTD_Z(MAXQ),    QPOSECO(MAXQ),
     &                  QTRK_X(MAXQ),     QTRK_Y(MAXQ),
     &                  QTRK_Z(MAXQ),     QTRK_P(MAXQ),
     &                  QTRK_N(MAXQ),     QTRK_D(MAXQ),
     &                  QTRK_C(MAXQ),     QTRK_T(MAXQ),
     &                  QTRK_V(MAXQ),     QTRK_I(MAXQ),
     &                  QNR_GOOD

      INTEGER           QEDCNCAND, QNR_GOOD
      REAL              EMPZ_TOT,EEMC_TOT,EHAD_TOT,PT_TOT,
     &                  FCALE50,  FCALE70,  FCALE90, FCALEm,
     &                  FCALEE50, FCALEE70, FCALEE90,FCALEEm,
     &                  QEPROB,         QECELNR,
     &                  QCAL_EE,        QCAL_XP,     QCAL_YP,
     &                  QCAL_ZP,
     &                  QEInList,       QENotInList,
     &                  QSRTD_E,        QSRTD_X,
     &                  QSRTD_Y,
     &                  QSRTD_Z,        QPOSECO,
     &                  QTRK_X,         QTRK_Y,
     &                  QTRK_Z,         QTRK_P,
     &                  QTRK_N,         QTRK_D,
     &                  QTRK_C,         QTRK_T,
     &                  QTRK_V,         QTRK_I

      CHARACTER*(650)   QEDTUP_form(2)
      DATA              QEDTUP_form /
     &                  'QEDTUP',
     &'QEDCNCAND:I,EMPZ_TOT:R,EEMC_TOT:R,EHAD_TOT:R,PT_TOT:R,FCALE50:R,
     &FCALE70:R,FCALE90:R,FCALEm:R,FCALEE50:R,FCALEE70:R,FCALEE90:R,
     &FCALEEm,QEPROB(2):R,QECELNR(2):R,QCAL_EE(2):R
     &,QCAL_XP(2):R,QCAL_YP(2):R,QCAL_ZP(2):R,QEInList(2):R,
     &QENotInList(2):R,QSRTD_E(2):R,QSRTD_X(2):R,QSRTD_Y(2):R,
     &QSRTD_Z(2):R,QPOSECO(2):R,QTRK_X(2):R,QTRK_Y(2):R,QTRK_Z(2):R,
     &QTRK_P(2):R,QTRK_N(2):R,QTRK_D(2):R,QTRK_C(2):R,QTRK_T(2):R,
     &QTRK_V(2):R,QTRK_I(2):R,QNR_GOOD:I'/






