      INTEGER           TLT_CLAS, TLT_SUBT, DST_BITS, TRIGB
      COMMON / TRGTUP / TLT_CLAS, TLT_SUBT, DST_BITS, TRIGB(4)

c -------------------------
c FLT BITS
c -------------------------
      LOGICAL FLTF,FLTP
      REAL TrigDat
      COMMON / FLTTUP / TrigDat(15)

C ----------------
C --- SLT bits ---
C ----------------
      real              SLT_EMPZ
      COMMON / TRGTP2 / SLT_EMPZ

C ----------------
C --- TLT bits ---
C ----------------
      integer TLT

c      TLTSUB1(32),  TLTSUB3(32),  TLTSUB4(32),
c     &                  TLTSUB5(32),  TLTSUB6(32),  TLTSUB7(32),
c     &                  TLTSUB8(32),  TLTSUB9(32),  TLTSUB10(32),
c     &                  TLTSUB11(32), TLTSUB12(32), TLTSUB13(32)
c      LOGICAL           TLTSUB1,      TLTSUB3,      TLTSUB4,
c     &                  TLTSUB5,      TLTSUB6,      TLTSUB7,
c     &                  TLTSUB8,      TLTSUB9,      TLTSUB10,
c     &                  TLTSUB11,     TLTSUB12,     TLTSUB13

      COMMON / TLTTUP / TLT(15)

c -----------------
c Tracking
c -----------------
      REAL             VCT_XV,  VCT_YV,  VCT_ZV,  CHVC,
     &                 VCT_XVC, VCT_YVC, VCT_ZVC, CHVCC,
     &                 FCAL_VTX, FCAL_VTXE, FTD_VTX
      INTEGER          NVTRK, NTRK, NVTRKC, NTRKC, FTD_VSEG

      COMMON / TRKTUP /VCT_XVC,VCT_YVC,VCT_ZVC,NVTRKC,NTRKC,CHVCC,
     &                 FCAL_VTX, FCAL_VTXE, FTD_VTX, FTD_VSEG

c ----------------------
c Calorimeter
c ----------------------
      REAL              FEMC_EN, BEMC_EN, REMC_EN,
     &                  FHAC_EN, BHAC_EN, RHAC_EN,
     &                  ET_TOT, calempz,  ETAMAX,
     &                  FCAL_TM, RCAL_TM,
     &                  BCAL_TM, 
     &                  GLOB_TM, cal_tot, PZ_TOT2,
     &                  hadempz,CPX, CPY, CPZ, CET

      COMMON / CALTUP / FEMC_EN, BEMC_EN, REMC_EN,
     &                  FHAC_EN, BHAC_EN, RHAC_EN,
     &                  FCAL_TM, RCAL_TM,
     &                  BCAL_TM, 
     &                  GLOB_TM, CPX, CPY, CPZ, CET, ETAMAX

      COMMON / CALTP1 / Fesum, Resum, Gesum, Uesum, Desum,
     &                  NFpmt, NRpmt, NGpmt, NUpmt, NDpmt,
     &                  Gerr,  Ferr,  Rerr,  FRerr, UDerr

      REAL              Fesum, Resum, Gesum, Uesum, Desum,
     &                  Gerr,  Ferr,  Rerr,  FRerr, UDerr

      INTEGER           NFpmt, NRpmt, NGpmt, NUpmt, NDpmt

c ------------------------
c Electron
c ------------------------
      INTEGER    MAXE
      PARAMETER (MAXE = 3)

      INTEGER           NE, IEM, IS1, IS4, CANDID, ICAND

      REAL              X_E,           Q2_E,       Y_E,
     &     EPROB,
     &     CAL_EE,        CAL_XP,     CAL_YP,
     &     CAL_ZP,
     &     EInList,       ENotInList,
     &     EMCenrgy,      HACenrgy,   Relec,
     &     SRTD_E,        SRTD_X,
     &     SRTD_Y,
     &     SRTD_Z,        POSECO,
     &     PRSE1,         PRSE2,
     &     PRSE3,
     &     TRK_X,         TRK_Y,
     &     TRK_Z,         TRK_P,
     &     TRK_N,         TRK_D,
     &     TRK_C,         TRK_T,
     &     TRK_V,         TRK_I,
     &     HES_X,         HES_Y,
     &     HES_E,         HES_R,
     &     HES_F,
     &     PXEL,          PYEL,
     &     PZEL,          ENEL,
     &     ECOR1,ECOR2, SRTD_ECOR
 
      INTEGER         ELCELLS, CellList, NrEcells,ecelnr

      COMMON / ELETUP / EPROB,
     &                  ECELNR,
     &                  CAL_EE,    CAL_XP, CAL_YP,
     &                  CAL_ZP,
     &                  EInList,   ENotInList,
     &                  EMCenrgy,  HACenrgy, Relec,
     &                  SRTD_E,    SRTD_X,
     &                  SRTD_Y,
     &                  SRTD_Z,    POSECO,
     &                  PRSE1,     PRSE2,
     &                  PRSE3,
     &                  TRK_X,     TRK_Y,
     &                  TRK_Z,     TRK_P,
     &                  TRK_N,     TRK_D,
     &                  TRK_C,     TRK_T,
     &                  TRK_V,     TRK_I,
     &                  HES_X,     HES_Y,
     &                  HES_E,     HES_R,
     &                  HES_F,
     &                  SRTD_ECOR, CANDID

      COMMON /ELETMP/ ELCELLS(50,MAXE),CellList(56),NrEcells

      INTEGER           NE1, Ncells
      
      REAL      ZPxBSp, ZPyBSp, ZPzBSp, ZPxCal, ZPyCal, ZPzCal, 
     &ZEminPz,ZufoCal,ZufoBSp,Zgamma,ZufoE,ZufoPx,ZufoPy,ZufoPz,
     &zuet,zupt

      COMMON /ZUFOTUP/ ZPxBSp, ZPyBSp, ZPzBSp,ZufoBSp,Zgamma,
     &     ZEminpz,zuet,zupt,
     &     ZPxCal(3,2), ZPyCal(3,2), ZPzCal(3,2), ZufoCal(3,2),
     &     ZufoE,ZufoPx,ZufoPy,ZufoPz

c ------------------------------------
c CorandCut
c ------------------------------------
       REAL empzcor, ptcor, etcor
       REAL y_jbcor, q2_jbcor, x_jbcor

       COMMON/HCOTUP/empzcor, ptcor, etcor, y_jbcor, q2_jbcor, 
     &      x_jbcor

c ------------------------------------
c LUMI
c ------------------------------------
      REAL POSX_LG, POSY_LG, ENE_LG,
     &     POSX_LE, POSY_LE, ENE_LE,
     &     TAG8_X,  TAG8_E,  TAG8_P,
     &     TAG8_V,
     &     TEMPLUMG,TEMPLUME,LUMI_XGTMP,LUMI_YGTMP,LUMI_XETMP, 
     &     LUMI_YETMP

      INTEGER   TAG44_E, ENE44M

      COMMON / LUMTUP / TEMPLUMG,TEMPLUME,ENE44M,TAG44_E,
     &     ENE_LE,POSX_LE,POSY_LE,
     &     ENE_LG,POSX_LG,POSY_LG

c --------------------------------------------
c General
c --------------------------------------------
      INTEGER           RUN_NUM, EV_NO, YEAR, FLT_NO,sparks
      REAL              EE, EP
      LOGICAL           montecarlo

      COMMON / GENTUP / RUN_NUM, EV_NO, YEAR, EE, EP,
     &     montecarlo,sparks
c ---------------------------------------------
c MCTRUE 
c ---------------------------------------------
      REAL              MCELE,    ET_MC,   YDA_MC,   XDA_MC,  Q2DA_MC,
     &                  MC_XV,    MC_YV,    MC_ZV,    SIGMA,    X_TRU,
     &                  Q2_TRU,   Y_TRU,           HCL_ENEG, HCL_PXEG,
     &                  HCL_PYEG, HCL_PZEG,EVT_WT, HCL_ENEE, HCL_PXEE,
     &                  HCL_PYEE, HCL_PZEE
      INTEGER           ICHN,genid
      COMMON / MR2TUP / MCELE,    ET_MC,   YDA_MC,   XDA_MC,  Q2DA_MC,
     &                  MC_XV,    MC_YV,    MC_ZV,    SIGMA,    X_TRU,
     &                  Q2_TRU,   Y_TRU,     ICHN, HCL_ENEG, HCL_PXEG,
     &                  HCL_PYEG, HCL_PZEG,EVT_WT, HCL_ENEE, HCL_PXEE,
     &                  HCL_PYEE, HCL_PZEE

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      INTEGER           EVTwant, EXOwant, XSECwant,
     &                  QEDC_NO, MUON_NO, C5_NO, VWALL_NO
 
      COMMON / BGDTUP / EVTwant, EXOwant, XSECwant,
     &                  QEDC_NO, MUON_NO, C5_NO, VWALL_NO
 

      CHARACTER*(500)   flt_form(2)
      DATA         flt_form /'FLTTUP',
     &'TrigDat(15):R' /
 
      CHARACTER*(500)   tlt_form(2)
      DATA         tlt_form /'TLTTUP',
     &'TLT(15):I' /

      CHARACTER*(500)    trk_form(2)
      DATA             trk_form /'TRKTUP',
     &'VCT_XVC:R,VCT_YVC:R,VCT_ZVC:R,NVTRKC:I,NTRKC:I,CHVCC:R,
     &FCAL_VTX:R,FCAL_VTXE:R,FTD_VTX:R,FTD_VSEG:I'/

      CHARACTER*(500)  cal_form(2)
      DATA              cal_form /'CALTUP',
     &'FEMC_EN:R,BEMC_EN:R,REMC_EN:R,FHAC_EN:R,BHAC_EN:R,
     &RHAC_EN:R,FCAL_TM:R,RCAL_TM:R,BCAL_TM:R,GLOB_TM:R,
     &CPX:R,CPY:R,CPZ:R,CET:R,ETAMAX:R' /

      CHARACTER*(500)   elec_form(2)
      DATA              elec_form /'ELETUP',
     &'EPROB:R,ECELNR:R,CAL_EE:R,CAL_XP:R,
     &CAL_YP:R,CAL_ZP:R,EInList:R,ENotInList:R,
     &EMCenrgy:R,HACenrgy:R,Relec:R,SRTD_E:R,SRTD_X
     &:R,SRTD_Y:R,SRTD_Z:R,POSECO:R,PRSE1:R,PRSE2:R
     &,PRSE3:R,TRK_X:R,TRK_Y:R,TRK_Z:R,TRK_P:R,TRK_
     &N:R,TRK_D:R,TRK_C:R,TRK_T:R,TRK_V:R,TRK_I
     &:R,HES_X:R,HES_Y:R,HES_E:R,HES_R:R,HES_F:R,
     &SRTD_ECOR:R,CANDID:I'/

      CHARACTER*(500)   zufos_form(2)
      DATA              zufos_form /'ZUFOTUP',
     &	   'ZPxBSp:R,ZPyBSp:R,ZPzBSp:R,ZufoBSp:R,Zgamma:R,
     &ZEminPz:R,Zuet:R,Zupt:R,ZPxCal(3,2):R,ZPyCal(3,2):R,
     &ZPzCal(3,2):R,ZufoCal(3,2):R,ZufoE:R,ZufoPx:R,ZufoPy:R,
     &ZufoPz:R'/
  
      CHARACTER*(500)    lumi_form(2) 
      DATA              lumi_form /'LUMTUP',
     &'TEMPLUMG:R,TEMPLUME:R,ENE44M:I,TAG44_E:I,
     &ENE_LE:R,POSX_LE:R,POSY_LE:R,
     &ENE_LG:R,POSX_LG:R,POSY_LG:R'/

      CHARACTER*(500)   gen_form(2)
      DATA              gen_form /'GENTUP',
     &    'RUN_NUM:I,EV_NO:I,YEAR:I,EE:R,EP:R,
     &    montecarlo:L,sparks:I' /

      CHARACTER*(500)	mctrue_form(2)
      DATA               mctrue_form /'MR2TUP',
     &'MCELE:R,ET_MC:R,YDA_MC:R,XDA_MC:R,Q2DA_MC:R,MC_XV:R,MC_YV:R,
     &MC_ZV:R,SIGMA:R,X_TRU:R,Q2_TRU:R,Y_TRU:R,ICHN:I,HCL_ENEG:R,
     &HCL_PXEG:R,HCL_PYEG:R,HCL_PZEG:R,EVT_WT:R,HCL_ENEE:R,
     &HCL_PXEE:R,HCL_PYEE:R,HCL_PZEE:R'/

      CHARACTER*(500)   BGDTUP_form(2)
      DATA              BGDTUP_form /
     &                  'BGDTUP',
     &'EVTwant:I,EXOwant:I,XSECwant:I,
     &QEDC_NO:I,MUON_NO:I,C5_NO:I,VWALL_NO:I'/












