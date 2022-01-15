      Subroutine Rewt_Vtx_Isr( Zvtx_true, Icheck,
     +     VtxWeight, Ierr ) 
      Implicit None
C     -------------------------------------------
C     INPUT:  Zvtx_true: Generated Z vertex
C             Icheck   : 0     Default
C                        +/- 1 Up/Down negative satellite
C                        +/- 2 Up/Down postiveive satellite
C     OUTPUT: VtxWeight
C
C     Version v1.0
C     N.Tuning 15/05/00
C     ------------------------
C
C     Error flag: If .ne. 0 Then vtxweight=1.
C
C     1: |Zvtx|>100 cm ( ==> wgt=1 )
C     2: Reweight from unknown mc vertex distribution.
C     3: Not reweight to 1 (high Q2 = 25cm) or 2 (low Q2 = 12x6)
C     5: 2+3
C     ----------------------------------------------------------
      Real    Zvtx_true, VtxWeight
      Real    VtxWeight2
      Integer Icheck
      Integer Ierr, Ierr2
      Logical Debug
      Common  /VtxNiels/ Debug
      Logical First
      Data    First/ .TRUE./
C     --------------------------------------------------------
c      Debug = .true.
	Debug = .false.
      If (First) Then
         First = .FALSE.
         Write(*,*) '=========================================='
         Write(*,*) ' This is Vertex Reweighting (ISR)         '
         Write(*,*) '             Version v1.0                 '
         Write(*,*) '                                          '
         Write(*,*) ' Date 15/05/2000 (N.Tuning)               '
         Write(*,*) ' Only to be used for 1996 ISR analysis    '
         Write(*,*) '=========================================='
      EndIf         
      If (Debug) Then
         Write(*,*) 
         Write(*,*) '------- Entering Vertex reweighting -----------'
         Write(*,*) '(Input:',Zvtx_true , Icheck,')'
      EndIf

      Ierr       = 0
      VtxWeight2 = 1.
      VtxWeight  = 1.

      If ( ABS(Zvtx_true).GE.100.) Then
         Ierr = 1
C         Write(*,*) 'WARNING Rewt_Vtx: |Zvtx|>100 cm:',Zvtx_true
         VtxWeight = 1.
         Return
      Endif

C     -----------------------------------
C     - Reweight MC to the DATA:
C     -----------------------------------
      Call Rewt_mc_to_data( Zvtx_true, Icheck,
     +     VtxWeight2, Ierr2) 


      VtxWeight = VtxWeight2
      Ierr      = Ierr+Ierr2

      If (Ierr .ne. 0) Then
         If (VtxWeight .ne.1) 
     +        Write(*,*) 'Warning error, but vtxwgt<>1:',ierr,VtxWeight
         VtxWeight = 1.
      EndIf

      If (Debug) Then
         Write(*,*)
         Write(*,*) '======== TOTAL VTX WEIGHT ============'
         Write(*,*) '======== ',VtxWeight,'    ============'
         Write(*,*) '======== ',Icheck,'               ============'
         Write(*,*) '======================================'
         Write(*,*)
      EndIf

      End

C     ****************************************************************
C     ****************************************************************


C     ==================================================================
      Subroutine Rewt_mc_to_data( Zvtx_true, Icheck,
     +     VtxWeight2, Ierr2) 
C     ==================================================================
c     Input: Zvtx    : Generated Zvtx
c     Output:VtxWeight2 : Weight
c            Ierr2      : Error flag
c     -------------------------------------------------------------
      Implicit None
      Real    Zvtx_true, VtxWeight2, X
      Integer Icheck
      Integer Ierr2
      Integer I
      Real    parmc(5),widthmc,meanmc,prbumc(2)
      Real    parda(5),widthda,meanda,prbuda(2)
C     ------------------------------------------------------
      Logical Debug
      Common  /VtxNiels/ Debug
      Real    IntF
      Data    IntF/2.5066283/        ! sqrt(2.*pi)
C     ------------------------------------------------------
      Real gauda(5),bunchda(5)
      Real gaumc(5),bunchmc(5),parmcv(5)
      Real DataInt, McInt, Normalize
      Real Vtxda, Vtxmc
C---
      Data parda/
     +     5.30501,
     +     1.67537,
     +     120.707,
     +     1.61893,
     +     2.9837 /
      Data widthda / 11.0539/
      Data meanda  / -3.60363/
      Data prbuda  / 67.2566,43.1767/
      Data parmc/
     +     11.4883,
     +     12.3699,
     +     409.947,
     +     4.90515,
     +     12.1458 /
      Data widthmc / 10.911/
      Data meanmc  / -3.19909/
      Data prbumc  / 69.0216,35.0813/
C     -------------------------------------------
      VtxWeight2 = 1.
      Ierr2      = 0
	debug = .false.
C     =============Initialize================
      Call Vzero(gauda,5)
      Call Vzero(gaumc,5)

      bunchda(1) = -1.*prbuda(1)
      bunchda(2) = -1.*prbuda(2)
      bunchda(3) = 0.
      bunchda(4) = prbuda(2)
      bunchda(5) = prbuda(1)

      bunchmc(1) = -1.*prbumc(1)
      bunchmc(2) = -1.*prbumc(2)
      bunchmc(3) = 0.
      bunchmc(4) = prbumc(2)
      bunchmc(5) = prbumc(1)

      parmcv(1)=parmc(1)
      parmcv(2)=parmc(2)
      parmcv(3)=parmc(3)
      parmcv(4)=parmc(4)
      parmcv(5)=parmc(5)
C     ================================

C----------------------------
C-- Vary Height of Satelites:
C----------------------------
      If (Icheck.eq.0) Then
C     DO NOTHING: DEFAULT
      ElseIf (Icheck.eq.-1) Then
         parmcv(1) = 0.70*parmc(1)
         parmcv(2) = 0.85*parmc(2)
      ElseIf (Icheck.eq. 1) Then
         parmcv(1) = 1.30*parmc(1)
         parmcv(2) = 1.15*parmc(2)
      ElseIf (Icheck.eq.-2) Then
         parmcv(5) = 0.70*parmc(5)
         parmcv(4) = 0.85*parmc(4)
      ElseIf (Icheck.eq. 2) Then
         parmcv(5) = 1.30*parmc(5)
         parmcv(4) = 1.15*parmc(4)
      EndIf


C     =============
      X = Zvtx_true
C     =============


C---------------------------------------------------------------
C-- Get Height of Vertex Distribution at X (sum of 5 gaussians):
C---------------------------------------------------------------
      Do i=1,5
         gauda(i) = parda(i)*exp(-0.5*
     +        ( (X-(meanda+bunchda(i)))/widthda) **2 )
         gaumc(i) = parmcv(i)*exp(-0.5*
     +        ( (X-(meanmc+bunchmc(i)))/widthmc) **2)
      Enddo
      Vtxda = gauda(1)+gauda(2)+gauda(3)+gauda(4)+gauda(5)
      Vtxmc = gaumc(1)+gaumc(2)+gaumc(3)+gaumc(4)+gaumc(5)


C--------------------------------------------------
C-- Normalize to integral of Vertex Distributions:
C--------------------------------------------------
      DataInt = ( widthda* (
     +     parda(1)+parda(2)+parda(3)+
     +     parda(4)+parda(5) 
     +     ) )
      McInt = ( widthmc* (
     +     parmcv(1)+parmcv(2)+parmcv(3)+
     +     parmcv(4)+parmcv(5) 
     +     ) )

      If (DataInt.gt.0.) Then
         Normalize = McInt/DataInt
      Else
         Normalize = 1.
         Write(*,*) 'rewt_vtx_isr ERROR DataInt:',DataInt
      EndIf

C---------------------------------
C-- And finally the Vertex Weight:
C---------------------------------
C     =================================
      If (VtxMc.ge.0.) Then
         Vtxweight2 = Normalize*Vtxda/Vtxmc
      Else
         Write(*,*) 'rewt_vtx_isr ERROR Vtxmc',Vtxmc
         VtxWeight2 = 1.
      EndIf
C     =================================

      If (Debug) Then
         Write(*,*) '=== Weight DATA to MC with 5 Gaussians === '
         Write(*,*) 'Vertex weight at Z vertex =',Zvtx_true,' cm'
         Write(*,*) '----------------------------------------------'
         Write(*,*) 'Integral of vtx distr. (DATA)   :',DataInt*IntF
         Write(*,*) 'Integral of vtx distr. (MC)     :',McInt*IntF
         Write(*,*) 'Norm.factor (MC/DATA of integrals):',Normalize
         Write(*,*) 'Vtx height at Z vertex (DATA)   :',Vtxda
         Write(*,*) 'Vtx height at Z vertex (MC)     :',Vtxmc
         Write(*,*) 'Ratio of heights (DATA/MC)      :',Vtxda/Vtxmc
         Write(*,*) '----------------------------------------------'
         Write(*,*) 'Vertex Weight :',Vtxweight2
         Write(*,*) 'Vary Satelites:',Icheck
         Write(*,*)
      EndIf

      End













