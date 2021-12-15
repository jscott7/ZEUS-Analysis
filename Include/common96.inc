      COMMON / TRGTUP / TLT_CLAS, TLT_SUBT, DST_BITS, TRIGB(4)

      INTEGER           TLT_CLAS, TLT_SUBT, DST_BITS, TRIGB

c -------------------------
c FLT BITS
c -------------------------
      COMMON / TRGTP1 / FLTF(64), FLTP(64), TrigDat(15)
      LOGICAL FLTF,FLTP
      REAL TrigDat

C ----------------
C --- SLT bits ---
C ----------------
      COMMON / TRGTP2 / SLT_EMPZ

      real              SLT_EMPZ

C ----------------
C --- TLT bits ---
C ----------------
      COMMON / TRGTP3 / TLT(15)
      
      integer TLT

c      TLTSUB1(32),  TLTSUB3(32),  TLTSUB4(32),
c     &                  TLTSUB5(32),  TLTSUB6(32),  TLTSUB7(32),
c     &                  TLTSUB8(32),  TLTSUB9(32),  TLTSUB10(32),
c     &                  TLTSUB11(32), TLTSUB12(32), TLTSUB13(32)
c      LOGICAL           TLTSUB1,      TLTSUB3,      TLTSUB4,
c     &                  TLTSUB5,      TLTSUB6,      TLTSUB7,
c     &                  TLTSUB8,      TLTSUB9,      TLTSUB10,
c     &                  TLTSUB11,     TLTSUB12,     TLTSUB13

c -----------------
c Tracking
c -----------------
      COMMON / TRKTUP /VCT_XV,  VCT_YV,  VCT_ZV,  NVTRK,  NTRK, CHVC,
     &                 VCT_XVC, VCT_YVC, VCT_ZVC, NVTRKC, NTRKC,CHVCC,
     &                 FCAL_VTX, FCAL_VTXE, FTD_VTX, FTD_VSEG

      REAL             VCT_XV,  VCT_YV,  VCT_ZV,  CHVC,
     &                 VCT_XVC, VCT_YVC, VCT_ZVC, CHVCC,
     &                 FCAL_VTX, FCAL_VTXE, FTD_VTX
      INTEGER          NVTRK, NTRK, NVTRKC, NTRKC, FTD_VSEG

c ----------------------
c Calorimeter
c ----------------------
      COMMON / CALTUP / FEMC_EN, BEMC_EN, REMC_EN,
     &                  FHAC_EN, BHAC_EN, RHAC_EN,
     &                  FCAL_TM, RCAL_TM,
     &                  BCAL_TM, 
     &                  GLOB_TM, CPX, CPY, CPZ, CET, ETAMAX

      REAL              FEMC_EN, BEMC_EN, REMC_EN,
     &                  FHAC_EN, BHAC_EN, RHAC_EN,
     &                  ET_TOT, calempz,  ETAMAX,
     &                  FCAL_TM, RCAL_TM,
     &                  BCAL_TM, 
     &                  GLOB_TM, cal_tot, PZ_TOT2,
     &                  hadempz,CPX, CPY, CPZ, CET

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
      COMMON / ELETUP / X_E,             Q2_E,       Y_E,
     &                  NE,  IS1,  IS4, IEM, ICAND,
     &                  EPROB,
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
     &                  PXEL(3),      PYEL(3),
     &                  PZEL(3),      ENEL(3),
     &        ECOR1(3),  ECOR2(3), SRTD_ECOR, CANDID


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

      COMMON /ELETMP/ ELCELLS(50,MAXE),CellList(56),NrEcells
      INTEGER         ELCELLS, CellList, NrEcells,ecelnr

      COMMON /ELETUP1 / NE1,
     &     ZPxBSp, ZPyBSp, ZPzBSp,ZufoBSp,
     &     ZPxCal(3,2), ZPyCal(3,2), ZPzCal(3,2), ZufoCal(3,2),
     &     ZEminPz, Zgamma,
     &     ZufoE,ZufoPx,ZufoPy,ZufoPz,
     &     Zupt,zuet, Ncells

       INTEGER           NE1, Ncells
      REAL      ZPxBSp, ZPyBSp, ZPzBSp, ZPxCal, ZPyCal, ZPzCal, 
     &ZEminPz,ZufoCal,ZufoBSp,Zgamma,ZufoE,ZufoPx,ZufoPy,ZufoPz,
     &zuet,zupt

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
      COMMON / LUMTUP / POSX_LG, POSY_LG, ENE_LG,
     &     POSX_LE, POSY_LE, ENE_LE,
     &     TAG8_X,  TAG8_E,  TAG8_P,
     &     TAG8_V,
     &     TAG44_E, 
     &     TEMPLUMG,TEMPLUME,LUMI_XGTMP,
     &     LUMI_YGTMP,LUMI_XETMP,LUMI_YETMP,
     &     ENE44M

      REAL POSX_LG, POSY_LG, ENE_LG,
     &     POSX_LE, POSY_LE, ENE_LE,
     &     TAG8_X,  TAG8_E,  TAG8_P,
     &     TAG8_V,
     &     TEMPLUMG,TEMPLUME,LUMI_XGTMP,LUMI_YGTMP,LUMI_XETMP, 
     &     LUMI_YETMP

      INTEGER   TAG44_E, ENE44M

c --------------------------------------------
c General
c --------------------------------------------
      COMMON / GENTUP / RUN_NUM, EV_NO, YEAR, FLT_NO, EE, EP,
     &     montecarlo,sparks

      INTEGER           RUN_NUM, EV_NO, YEAR, FLT_NO,sparks
      REAL              EE, EP
      LOGICAL           montecarlo

c ---------------------------------------------
c MCTRUE 
c ---------------------------------------------
      COMMON / MR2TUP / MCELE,    ET_MC,   YDA_MC,   XDA_MC,  Q2DA_MC,
     &                  MC_XV,    MC_YV,    MC_ZV,    SIGMA,    X_TRU,
     &                  Q2_TRU,   Y_TRU,     ICHN, HCL_ENEG, HCL_PXEG,
     &                  HCL_PYEG, HCL_PZEG,EVT_WT, HCL_ENEE, HCL_PXEE,
     &                  HCL_PYEE, HCL_PZEE,genid
      REAL              MCELE,    ET_MC,   YDA_MC,   XDA_MC,  Q2DA_MC,
     &                  MC_XV,    MC_YV,    MC_ZV,    SIGMA,    X_TRU,
     &                  Q2_TRU,   Y_TRU,           HCL_ENEG, HCL_PXEG,
     &                  HCL_PYEG, HCL_PZEG,EVT_WT, HCL_ENEE, HCL_PXEE,
     &                  HCL_PYEE, HCL_PZEE
      INTEGER           ICHN,genid

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      INTEGER           EVTwant, EXOwant, XSECwant,
     &                  QEDC_NO, MUON_NO, C5_NO, VWALL_NO
 
      COMMON / BGDTUP / EVTwant, EXOwant, XSECwant,
     &                  QEDC_NO, MUON_NO, C5_NO, VWALL_NO
 

      CHARACTER*(500)   flt_form(2)
      DATA         flt_form /'FLT',
     &'TrigDat(15):R' /
 
      CHARACTER*(500)   tlt_form(2)
      DATA         tlt_form /'TLT',
     &'TLT(15):I' /

      CHARACTER*(500)    trk_form(2)
      DATA             trk_form /'TRK',
     &'VCT_XVC:R,VCT_YVC:R,VCT_ZVC:R,NVTRKC:I,NTRKC:I,CHVCC:R,
     &FCAL_VTX:R,FCAL_VTXE:R,FTD_VTX:R,FTD_VSEG:I'/

      CHARACTER*(500)  cal_form(2)
      DATA              cal_form /'CAL',
     &'FEMC_EN:R,BEMC_EN:R,REMC_EN:R,FHAC_EN:R,BHAC_EN:R,
     &RHAC_EN:R,FCAL_TM:R,RCAL_TM:R,BCAL_TM:R,GLOB_TM:R,
     &CPX:R,CPY:R,CPZ:R,CET:R,ETAMAX:R' /

      CHARACTER*(500)   elec_form(2)
      DATA              elec_form /'ELEC',
     &'EPROB:R,ECELNR:R,CAL_EE:R,CAL_XP:R,
     &CAL_YP:R,CAL_ZP:R,EInList:R,ENotInList:R,
     &EMCenrgy:R,HACenrgy:R,Relec:R,SRTD_E:R,SRTD_X
     &:R,SRTD_Y:R,SRTD_Z:R,POSECO:R,PRSE1:R,PRSE2:R
     &,PRSE3:R,TRK_X:R,TRK_Y:R,TRK_Z:R,TRK_P:R,TRK_
     &N:R,TRK_D:R,TRK_C:R,TRK_T:R,TRK_V:R,TRK_I
     &:R,HES_X:R,HES_Y:R,HES_E:R,HES_R:R,HES_F:R,
     &SRTD_ECOR:R,CANDID:I'/

      CHARACTER*(500)   zufos1_form(2)
      DATA              zufos1_form /'ZUFOS1',
     &	   'ZPxBSp:R,ZPyBSp:R,ZPzBSp:R,ZufoBSp:R'/

      CHARACTER*(500) zufos2_form(2)
      DATA           zufos2_form /'ZUFOS2','Zgamma:R'/

      CHARACTER*(500) zufos3_form(2)
      DATA           zufos3_form /'ZUFOS3','ZEminPz:R,
     &Zuet:R,Zupt:R'/

      CHARACTER*(500) zufos4_form(2)
      DATA           zufos4_form /'ZUFOS4',
     &     'ZPxCal(3,2):R,ZPyCal(3,2):R,ZPzCal(3,2):R, 
     &     ZufoCal(3,2):R'/

      CHARACTER*(500) zufos5_form(2)
      DATA           zufos5_form /'ZUFOS5',
     &     'ZufoE:R,ZufoPx:R,ZufoPy:R,ZufoPz:R'/

      CHARACTER*(500)    temp_form(2) 
      DATA              temp_form /'TEMP',
     &'TEMPLUMG:R,TEMPLUME:R'/

      CHARACTER*(500)    tag1_form(2)
      DATA 		tag1_form /'TAG1','ENE44M:I'/

      CHARACTER*(500)    tag2_form(2)
      DATA              tag2_form /'TAG2','TAG44_E:I'/

      CHARACTER*(500)    lumi_form(2) 
      DATA              lumi_form /'LUMI1',
     &'ENE_LE:R,POSX_LE:R,POSY_LE:R'/

      CHARACTER*(500)    lumi2_form(2)
      DATA              lumi2_form /'LUMI2',
     &'ENE_LG:R,POSX_LG:R,POSY_LG:R'/

      CHARACTER*(500)   gen_form(2)
      DATA              gen_form /'GEN',
     &    'RUN_NUM:I,EV_NO:I,YEAR:I,FLT_NO:I,EE:R,EP:R,
     &    montecarlo:L,sparks:I' /

      CHARACTER*(500)	mctrue_form(2)
      DATA               mctrue_form /'MCTRUE',
     &'MCELE:R,ET_MC:R,YDA_MC:R,XDA_MC:R,Q2DA_MC:R,MC_XV:R,MC_YV:R,
     &MC_ZV:R,SIGMA:R,X_TRU:R,Q2_TRU:R,Y_TRU:R,ICHN:I,HCL_ENEG:R,
     &HCL_PXEG:R,HCL_PYEG:R,HCL_PZEG:R,EVT_WT:R,HCL_ENEE:R,
     &HCL_PXEE:R,HCL_PYEE:R,HCL_PZEE:R'/

      CHARACTER*(500)   BGDTUP_form(2)
      DATA              BGDTUP_form /
     &                  'BGDTUP',
     &'EVTwant:I,EXOwant:I,XSECwant:I,
     &QEDC_NO:I,MUON_NO:I,C5_NO:I,VWALL_NO:I'/












