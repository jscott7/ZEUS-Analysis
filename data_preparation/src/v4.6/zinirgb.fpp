      Subroutine z_ini_RecGB(MoreTr, SecTr, EneCor, BSpCor)
      Implicit NONE
*-------------------------------------------------------------------------------
*
*   Called by USER in ZUINIT to initialize parameters for z_RecGB routine
*
*     INPUT:  NONE
*     OUTPUT: NONE
*
*   Author:   Gennady Briskin, Tel Aviv University
*   Date:     20-Dec-1996
*
*_______________________________________________________________________________
*
#include "zdrecgb.inc"
*
      Integer MoreTr, SecTr, EneCor, BSpCor
      Logical First
      DATA    First /.TRUE./
*
       AQ_options       = .FALSE.
       sec_vtx_yes      = .FALSE.
       EcorNT           = .FALSE.
       BcorNT           = .FALSE.

       If (MoreTr.eq.1) AQ_options       = .TRUE.
       If (SecTr.eq.1)  sec_vtx_yes      = .TRUE.
       If (EneCor.eq.1) EcorNT           = .TRUE.
       If (BSpCor.eq.1) BcorNT           = .TRUE.

       printflag        = .FALSE.
       zRecGB_Debug     = .FALSE. 
       vmode            = 0
       ptMin_CUT        = 0.1 
       matPrCut         = 0.02 !obsolete; matching program MATSCH

       If (EcorNt) Then  !!! corandcut values 
c          PrCutHad         = 0.1
c          PrCutElec        = 0.3
          PrCutHad         = 0.05
          PrCutElec        = 0.2
       Else
          PrCutHad         = 0.05
          PrCutElec        = 0.2
       Endif

       IF (AQ_options) THEN
          sl_CUT           = 4
          muon_eovp        = 0.25
          muon_Emax        = 5.0
          muon_Ptmax       = 30.
          EovPCut          = 1.0
          zNsig            = 1.2
          ptMax_CUT9       = 25.
          ptMax_CUT7       = 20.
       ELSE
          sl_CUT           = 3
          EovPCut          = 0.8
          zNsig            = 1.0
          ptMax_CUT9       = 15.
          ptMax_CUT7       = 15.
       ENDIF
       eCut_DCA         = 25.0  
       DCA_hCut_min     = 20.0 
C      EndIf
*
       If (First) Then
          First = .FALSE.
          Write(6,100)zRecGB_Debug,vmode,sl_CUT,ptMin_CUT,ptMax_CUT,
     &            PrCutElec,PrCutHad,matPrCut,EovPCut,
     &            zNsig,eCut_DCA,DCA_hCut_min,sec_vtx_yes
*
       EndIf

 100  Format(
     &10x,'*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'/
     &10x,'**                                          **'/
     &10x,'**     Parameters for z_RecGB routine       **'/
     &10x,'**                                          **'/
     &10x,'**      zRecGB_Debug = ',L8  ,'             **'/
     &10x,'**      matsch_vmode = ',I8  ,'             **'/
     &10x,'**      sl_CUT       = ',I8  ,'             **'/
     &10x,'**      ptMin_CUT    = ',F8.4,'             **'/
     &10x,'**      ptMax_CUT    = ',F8.4,'             **'/
     &10x,'**      PrCutElec    = ',F8.4,'             **'/
     &10x,'**      PrCutHad     = ',F8.4,'             **'/
     &10x,'**      matPrCut     = ',F8.4,'             **'/
     &10x,'**      EovPCut      = ',F8.4,'             **'/
     &10x,'**      zNsig        = ',F8.4,'             **'/
     &10x,'**      eCut_DCA     = ',F8.4,'             **'/
     &10x,'**      DCA_hCut_min = ',F8.4,'             **'/
     &10x,'**      sec_vtx_yes  = ',L8  ,'             **'/      
     &10x,'**                                          **'/
     &10x,'*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'/
     &      ) 
*
      End

