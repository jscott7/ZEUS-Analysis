
*==>  ZUFOs.....................................................
*...
*     Nzufos            <==> Number of ZUFOs
*     tufo(4,1:Nzufos)  <==> Type of ZUFO information
*     zufo(4,1:Nzufos)  <==> 4 momentum of ZUFO (1=px,2=py,3=pz,4=E)
*...
      Integer   zufos_Max
      Parameter(zufos_Max=250)
      Integer        Nzufos,tufo
      Real                                    zufo
      Common/zRec01/ Nzufos,tufo(4,zufos_Max),zufo(4,zufos_Max)

*==>  Electron Calorimeter Info.................................
      Integer         eNdis
      Real                  ePrdis,eEdis,eXdis,eYdis,eZdis
      Common/zRec02/  eNdis,ePrdis,eEdis,eXdis,eYdis,eZdis

*==>  Electro Track Info........................................
      Integer   neT_Max
      Parameter(neT_Max=10)      
      Integer         eTrMat, eTrNear, vcteID
      Real                                              eTrDCA
      Common/zRec03/  eTrMat, eTrNear, vcteID(neT_Max), eTrDCA

*==>  Hadronic Track Info.......................................
      Integer   nT_Max
      Parameter(nT_Max=100)      
      Integer         nT ,vcthID, qTr, swmTr, qltTr, mTrI, wTrI
      Real            pTr,  dpTr, thTr, phTr
      Real            Dist_IsTr
      Common/zRec04/ nT , vcthID(nT_Max), qTr(  nT_Max), swmTr(nT_Max), 
     &                      qltTr(nT_Max), mTrI( nT_Max),
     &                      pTr(  nT_Max), dpTr( nT_Max),
     &                      thTr( nT_Max), phTr( nT_Max),
     &                    wTrI(10,nT_Max)
     &                   ,Dist_IsTr( nT_Max)

*==>  Parameters................................................
      Logical  zRecGB_Debug
      Integer  vmode
      Integer  sl_CUT
      Real     PrCutHad             
      Real     PrCutElec
      Real     ptMin_CUT
      Real     ptMax_CUT
      Real     matPrCut
      Real     EovPCut
      Real     zNsig
      Real     eCut_DCA
      Real     DCA_hCut_min

      Common /zRec05/ zRecGB_Debug,vmode,sl_CUT,
     &                PrCutHad,PrCutElec,
     &                ptMin_CUT,ptMax_CUT,matPrCut,EovPCut,zNsig,
     &                eCut_DCA,DCA_hCut_min


*==> Additions; AQ 12/98

      Logical  AQ_options
      Real     ptMax_CUT7,ptMax_CUT9
      Real     muon_eovp
      Real     muon_Emax
      Real     muon_Ptmax

*-->  Switch added for including sec-vtz-tracks 7.10.98 UW
      logical  sec_vtx_yes, printflag

      Logical  EcorNT, BcorNT

      Common /zRec06/ AQ_options,ptMax_CUT7,ptMax_CUT9,
     &                muon_eovp,muon_emax,muon_ptmax,
     &                sec_vtx_yes, printflag, EcorNT, BcorNT
