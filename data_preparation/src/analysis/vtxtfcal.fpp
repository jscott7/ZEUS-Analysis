
C     =========================================
      Real Function vtxtfcal(NrGoodCells,VzErr)
C     =========================================
      Implicit None
C
      Integer NrGoodCells
      Real    VzErr
      Integer FNrCell,MaxFNrCell
      Parameter(MaxFNrCell=500)
C
C random generator variables
      INTEGER         NTOTIN,NTO2IN
      INTEGER         islnd,islnt,seed,igenno
      INTEGER         isl(40)
      REAL            rgenno
      COMMON /SLATE/  isl
C
      Logical First
      Save    First
      Data    First /.TRUE./
C
      Integer FCellIcel(MaxFNrCell),FCellIreg(MaxFNrCell),
     +        FCellImod(MaxFNrCell),FCellItow(MaxFNrCell),
     +        FCellId(MaxFNrCell)
      Real    FCellE(MaxFNrCell),FCellImb(MaxFNrCell),
     +        FCellT(MaxFNrCell),FCellTcor(MaxFNrCell),
     +        FCellIzol(MaxFNrCell),FEfcal,
     +        FCellWag(MaxFNrCell),FCellX(MaxFNrCell),
     +        FCellY(MaxFNrCell),FCellZ(MaxFNrCell),
     +        FCellVtx(MaxFNrCell),
     +        Ftint,FZint,FEfcalEMC
      Integer Fyear,Frun
      Common /VTXTFCALCOMMON/ FNrCell,FCellIcel,FCellIreg,
     +                        FCellImod,FCellItow,
     +                        FCellE,FCellImb,FCellWag,
     +                        FCellT,FCellTcor,FCellIzol,
     +                        FCellX,FCellY,FCellZ,FCellVtx,
     +                        Fyear,Frun,FCellId,FEfcal,
     +                        Ftint,FZint,FEfcalEMC
C
      Integer Ngood
      Logical MC
      Integer i,j,k
      Real emctab(-1:26,-1:101),
     +     hac1tab(-1:26,-1:26),
     +     hac2tab(-1:26,-1:26)
      Common /IZOLTABLETFV/ emctab,hac1tab,hac2tab
      Real EconeIzol,TfcalCorr,TimeWag,GetVz,DVz,vtxEfcalcorr
C
      Real MCVtx
      Real RmsVtz,RmsVtz1,RmsVtzMax,Vtz,WagInit,Wag
      Integer imax
C not used:
C      Real MinRMS,MinFracWag,MinCellsNr
C      Parameter (MinRMS=7.5)
C      Parameter (MinFracWag=0.9)
C      Parameter (MinCellsNr=25)
C
      Integer NrCalCell,CellVec(5)
      Integer Ireg,Icel,Imod,Itow
      Equivalence (Ireg,CellVec(1))
      Equivalence (Icel,CellVec(3))
      Equivalence (Imod,CellVec(4))
      Equivalence (Itow,CellVec(5))
      Character*5 CellKind
      Real ClX,CLY,ClZ,ee1,ee2
      Logical Error
C
      Integer IndCT
      Save    IndCT
C
      Real TimeMaxMin
      Parameter (TimeMaxMin=12.)
C
#include "partap.inc"
#include "caltru.inc"
#include "ctime2.inc"
#include "zrevt.inc"
#include "fmcrun.inc"
#include "usgrun.inc"
#include "fmckin.inc"
#include "usgevt.inc"
#include "fmcvtx.inc"
C
C
      If(First) Then
         IndCT=GetInd(Caltru,'E')
C...Initialize random number generator
         call datime(islnd,islnt)
         seed = 25*(3600*isl(4)+60*isl(5)+isl(6))
         if (mod(seed,2).eq.0) seed=seed+1
         call rmarin(seed,NTOTIN,NTO2IN)
         Write(*,*) 'SEED =',seed
C
         First=.FALSE.
      EndIf
C
      VzErr = 1000000.0
C
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Filling of local tables
      If ( CouTab(FMCRUN).GT.0 .OR. CouTab(USGRUN).GT.0 .OR.
     +     CouTab(FMCKin).GT.0 .OR. CouTab(USGEVT).GT.0) Then
         MC = .TRUE.
         If(CouTab(FMCVtx).GT.0)Then
            Call FetTab(FMCVTX,ID,1)
            MCVtx = FMCVTX_R(3)
         Else
            MCVtx = 0.0
         EndIf
      Else
         MC = .FALSE.
      EndIf
C
      If(Coutab(ZREVT).le.0) Return
      Call Fettab(ZREVT,ID,1)
      Frun =  ZREVT_RunNr
      Fyear = ZREVT_time(1)/10000
      If(Coutab(CTIME2).eq.0) Then
          Fzint = 0.0
          Ftint = 0.0
      Else
          Call Fettab(CTIME2,ID,1)
          Fzint = CTIME2_zint
          Ftint = CTIME2_tint
      EndIf
      FEfcal  = 0.0
      FNrCell = 0
      FEfcalEMC = 0.0
      NrCalCell = Coutab(Caltru)
      Do i=NrCalCell,1,-1
         Call Fettab(Caltru,IndCT,i)
         Call CcDeco(CalTru_CellNr,CellKind,CellVec,Error)
         Call Cccxyz(CalTru_CellNr,ClX,ClY,ClZ,Error)
         If(Error) GoTo 1000
         If(Ireg.ne.0) GoTo 1000
         FEfcal = FEfcal + CalTru_E
         If(Icel.le.4) FEfcalEMC=FEfcalEMC+CalTru_E
         If((FNrCell.lt.MaxFNrCell).and.(CalTru_E.gt.0.001)) Then
            FNrCell = FNrCell + 1
            FCellIcel(FNrCell) = Icel
            FCellIreg(FNrCell) = Ireg
            FCellImod(FNrCell) = Imod
            FCellItow(FNrCell) = Itow
            FCellE(FNrCell)    = CalTru_E
            FCellImb(FNrCell)  = CalTru_Imbal
            ee1 = (CalTru_E - CalTru_Imbal)/2.
            ee2 = (CalTru_E + CalTru_Imbal)/2.
            If(ee1.gt.0.001) Then
                ee1 = 1/(0.40+1.4/(ee1**0.65))**2
            Else
                ee1 = 0.0
            EndIf
            If(ee2.gt.0.001) Then
                ee2 = 1/(0.40+1.4/(ee2**0.65))**2
            Else
                ee2 = 0.0
            EndIf
            FCellT(FNrCell)    = 1.E+9 *
     +                (CalTru_t(1)*ee1+CalTru_t(2)*ee2)/(ee1+ee2)
            FCellTcor(FNrCell) = 0.0
            FCellId(FNrCell)   = CalTru_ID
            FCellX(FNrCell)    = CLx
            FCellY(FNrCell)    = CLy
            FCellZ(FNrCell)    = CLz
         EndIf
C
 1000    Continue
      EndDo
C
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      If((FNrCell.lt.10).or.(FEfcal.lt.1.0)) Return
C
C
C     initialize table for izolation
C
      Call Vzero(emctab,2884)
      Call Vzero(hac1tab,784)
      Call Vzero(hac2tab,784)
C
      Do i=1,FNrCell
         If(FCellIreg(i).ne.0) GoTo 200
         If(FCellIcel(i).le.4) Then
            j = FCellImod(i)
            k = FCellICel(i)+FCellITow(i)*4
            emctab(j,k)=FCellE(i)
         ElseIf(FCellIcel(i).eq.6) Then
            j = FCellImod(i)
            k = FCellITow(i)
            hac1tab(j,k)=FCellE(i)
         ElseIf(FCellIcel(i).eq.7) Then
            j = FCellImod(i)
            k = FCellITow(i)
            hac2tab(j,k)=FCellE(i)
         EndIf
 200     Continue
      EndDo
C
      Ngood = 0
      Do i=1,FNrCell
         FCellIzol(i) = Econeizol( FCellImod(i),FCellItow(i)
     +                            ,FCellIcel(i),FCellIreg(i))
         FCellWag(i)  = TimeWag( FCellE(i),FCellIzol(i)
     +                          ,FCellImod(i),FCellItow(i)
     +                          ,FCellIcel(i),FCellIreg(i)
     +                          ,FCellImb(i))
         If(FCellWag(i).gt.0.0) Then
            FCellTcor(i) = TfcalCorr( FCellImod(i),FCellItow(i)
     +                            ,FCellIcel(i),FCellIreg(i)
     +                            ,Fyear,FCellE(i),FEfcal,Frun
     +                            ,FCellIzol(i),FCellT(i))
            If(abs(FCellTcor(i)).gt.TimeMaxMin) FCellWag(i)=0.0
         EndIf
C
         If(MC) Then
            If(FCellWag(i).gt.0.0) NGood  = NGood + 1
         Else
            If(FCellWag(i).gt.0.0) Then
               FCellVtx(i)  = getvz(FCellTcor(i),FCellX(FNrCell),
     +                              FCellY(FNrCell),FCellZ(FNrCell),
     +                              FZint,FTint)
               If(FCellVtx(i).eq.1000000.0) Then
                  FCellWag(i) = 0.0
               Else
                  NGood  = NGood + 1
               EndIf
            EndIf
C
         EndIf
      EndDo
C
      Wag = 0.0
      NrGoodCells = NGood
      If(NrGoodCells.eq.0) Return
      If(.NOT.MC) Then
         Vtz = 0.0
         Do i=1,FNrCell
            If(FCellWag(i).gt.0.0) Then
               Vtz = Vtz + FCellVtx(i)*FCellWag(i)
               Wag = Wag + FCellWag(i)
            EndIf
         EndDo
         Vtz = Vtz / Wag
         vtxtfcal = Vtz + vtxEfcalcorr(FEfcal,FEfcalEMC)
         VzErr    = DVz(NGood)
      Else
         VzErr    = DVz(NGood)
         Call rnorml(rgenno,1)
         vtxtfcal = rgenno * VzErr + MCVtx
      EndIf
C
      Return
      End



C     ===================================================
      Real Function TfcalCorr(imod,itow,icel,ireg,year,
     +                        Ecell,EFcal,run,eizol,time)
C     ===================================================
      Implicit None
C
C     Main Correction Function
C
C     Input  : imod/itow/icel/ireg - CAL element coordinates
C              year  - year ????
C              ecell - energy of cell
C              efcal - energy in FCAL
C              run   - run numer
C              eizol - izolation
C              time  - measured cell time
C
C     Output : TfcalCorr - corrected time (-1000.0 means error)
C
C
      Integer imod,itow,icel,ireg,year,run
      Real Ecell,EFcal,eizol,time
      Integer typ,izotyp
      Real TfEcellCorr,TfEfcalCorr,TfRunCorr,TfModCorr,TfEizolCorr
C
      TfcalCorr = -1000.0
      If(Ireg.ne.0) Return
      If(Icel.le.4) Then
         typ = 1
      ElseIf(Icel.eq.6) Then
         typ = 2
      ElseIf(Icel.eq.7) Then
         typ = 3
      Else
         Return
      EndIf
C
      If(typ.gt.1) Then
         If( (Imod.le.13).and.(Imod.ge.11).and.
     +       (Itow.le.13).and.(Itow.ge.11)      ) Then
             izotyp = 10*(typ-1)+1
         Else
             izotyp = 10*(typ-1)
         EndIf
      Else
         izotyp = 1
      EndIf
C
      TfcalCorr = ( time  +  TfEcellCorr(Ecell,typ)
     +                    +  TfEfcalCorr(Efcal,typ)
     +                    +  TfRunCorr(Run)
     +                    +  TfModCorr(imod,itow,icel,year)
     +                    +  TfEizolCorr(eizol,izotyp) ) / 0.862
C
      Return
      End



C     ================================
      Real Function TfEcellCorr(E,typ)
C     ================================
      Implicit None
C
C     Typ = 1 -> EMC
C
      Real E
      Integer typ
      Real A(0:9),EE
      Integer i
C
      Call Vzero(A,10)
C
      If(typ.eq.1) Then
         If(E.le.4.0) Then
            A(0)= -0.66271
            A(1) = 3.7359
            A(2) = -3.4849
            A(3) = 1.8276
            A(4) = -0.52554
            A(5) = 0.77199E-01
            A(6) = -0.45167E-02
         ElseIf(E.le.11.) Then
            A(0) = 1.1690
            A(1) = 0.97734E-01
            A(2) = -0.38792E-02
            A(3) = 0.54359E-04
         ElseIf(E.le.33.5) Then
             A(0) = 1.7352
             A(1) = 0.11192E-01
             A(2) = -0.11707E-03
         ElseIf(E.le.48.5) Then
             A(0) =  5.4193
             A(1) = -0.14873
             A(2) = 0.13599E-02
         ElseIf(E.le.120.0) Then
            A(0) = 1.0489
            A(1) = 0.80797E-02
         Else
            A(0) = 2.0
         EndIf
      ElseIf(Typ.eq.2) Then
         If(E.le.0.275) Then
            A(0) = -3.3
         ElseIf(E.le.3.2) Then
            A(0) = -5.9558
            A(1) = 13.619
            A(2) = -19.988
            A(3) = 17.142
            A(4) = -8.9309
            A(5) = 2.8577
            A(6) = -0.54756
            A(7) = 0.57561E-01
            A(8) = -0.25490E-02
         ElseIf(E.lt.9.) Then
            A(0)= -1.9517
            A(1)= 0.39484
            A(2)= -0.55602E-01
            A(3)= 0.42073E-02
            A(4) = -0.12424E-03
         ElseIf(E.lt.31.) Then
            A(0)= -1.0820
            A(1)= 0.54358E-01
            A(2)= -0.71507E-03
         ElseIf(E.lt.44.) Then
            A(0)= 2.0667
            A(1)= -0.97337E-01
            A(2)=  0.89439E-03
         ElseIf(E.lt.102.) Then
            A(0)= -0.99523
            A(1)= 0.11424E-01
         Else
            A(0)= 0.19629
         EndIf
      ElseIf(Typ.eq.3) Then
         If(E.lt.0.21) Then
            A(0)=-2.42
         ElseIf(E.lt.3.4) Then
            A(0)=-4.6685
            A(1)=13.741
            A(2)=-16.508
            A(3)=10.977
            A(4)= -4.2231
            A(5)= 0.93472
            A(6)=-0.11028
            A(7)= 0.53641E-02
         ElseIf(E.lt.10.5) Then
            A(0)=0.41140
            A(1)=0.14714
            A(2)=-0.54475E-02
         ElseIf(E.lt.20.95) Then
            A(0)=0.76501
            A(1)=0.69506E-01
            A(2)=-0.12666E-02
         Else
            A(0)=1.655
         EndIf
      EndIf
C
      TfEcellCorr=A(0)
      EE=1.0
      Do i=1,9
         EE=EE*E
         TfEcellCorr=TfEcellCorr+EE*A(i)
      EndDo
C
C      Write(*,*) ' Energy correction = ',TfEcellCorr
C
      Return
      End


C     =====================================
      Real Function TfEfcalCorr(efcal,tryb)
C     =====================================
      Implicit None
C
      Integer tryb
      Real efcal,x
C
      TfEfcalCorr = 0.0
      x = log10(efcal)
C
      If(tryb.eq.1) Then
        If(x.gt.0.5) Then
           Tfefcalcorr = 0.13584 + 0.17457 * x
     +                -0.14983 * x * x
        Else
           Tfefcalcorr = 0.187
        EndIf
      ElseIf(tryb.eq.2) Then
        If(x.gt.1.15) Then
           Tfefcalcorr = 1.0373  -0.51420 * x
        Else
           Tfefcalcorr = 0.45
        EndIf
      ElseIf(tryb.eq.3) Then
        If(x.gt.0.5) Then
           Tfefcalcorr = 1.1106  -0.55936 * x
        Else
           Tfefcalcorr = 0.8
        EndIf
      EndIf
C
C      Write(*,*) ' Efcal correction = ',Tfefcalcorr
C
      Return
      End


C     ============================================
      Real Function TfModCorr(imod,itow,icel,year)
C     ============================================
      Implicit None
C
      Integer imod,itow,icel,year
C
C
      Tfmodcorr = 0.0
C
      If(Year.ge.1995) Then
         If((IMod.eq.13).and.(ITow.eq.13).and.(icel.ne.3))Then
            TfModCorr = - 3.
            GoTo 333
         EndIf
         If((icel.eq.4).and.(ITow.eq.11).and.(IMod.eq.12))Then
            TfModCorr = 2.5
            GoTo 333
         EndIf
      ElseIf(Year.ge.1994) Then
         If((icel.eq.7).and.(ITow.eq.11).and.(IMod.eq.12))Then
            TfModCorr = + 2.5
            GoTo 333
         EndIf
      EndIf
C
C      Write(*,*) 'Modul correction = ',TfModCorr
C
 333  Continue
C
      Return
      End


C     ====================================
      Real Function TfEizolCorr(fizol,typ)
C     ====================================
      Implicit None
C
      Integer typ
      Real fizol
C
C     typ =  1  EMC
C           10  HAC1 1 ring
C           11  HAC1 - 1 ring
C           20  HAC2 1 ring
C           21  HAC2 - 1 ring
C
      TfEizolCorr = 0.0
C EMC
      If(typ.eq.1) Then
         If(fizol.gt.20.0) Then
            TfEizolcorr = -0.27154 + 0.18882E-01 * fizol
     +                    -0.31865E-03 * fizol**2
     +                    +0.17962E-05 * fizol**3
         ElseIf(fizol.gt.0.01) Then
            TfEizolcorr = -0.25990 + 0.12137E-01 * fizol
         EndIf
C HAC1 -ir
      ElseIf(typ.eq.10) Then
         TfEizolcorr = -3.8004
     +                 +5.2447*EXP(-4.2571*fizol**(-1.0028))
C HAC1 ir
      ElseIf(typ.eq.11) Then
         TfEizolcorr = -3.6339
     +                 +5.2378*EXP(-4.2604*fizol**(-0.81314))
C HAC2 -ir
      ElseIf(typ.eq.20) Then
         If(fizol.gt.2.3) Then
            TfEizolcorr = 1.1130-12.842*EXP(-0.81922*fizol**(0.40431))
         Else
            TfEizolcorr = -3.0
         EndIf
C HAC2 ir
      ElseIf(typ.eq.21) Then
         If(fizol.gt.45.) Then
            TfEizolcorr = -0.15492 + 0.10772E-01*fizol
         ElseIf(fizol.gt.20.) Then
            TfEizolcorr = -14.818
     +                  +3.8571*EXP(1.1870*fizol**(0.37557E-01))
         ElseIf(fizol.gt.2.) Then
            TfEizolcorr = -15.842
     +                  +3.7628*EXP(1.1942*fizol**(0.57753E-01))
         Else
            TfEizolcorr = -2.8
         EndIf
      EndIf
C
C      Write(*,*) 'Izolation correction  = ',TfEizolcorr
C
      Return
      End


C     ================================================
      Real Function Econeizol(iimod,iitow,iicel,iireg)
C     ================================================
      Implicit None
C
      Integer iimod,iitow,iicel,iireg
C
      Real emctab(-1:26,-1:101),
     +     hac1tab(-1:26,-1:26),
     +     hac2tab(-1:26,-1:26)
      Common /IZOLTABLETFV/ emctab,hac1tab,hac2tab
C
      Integer j,k,l,m
C
      econeizol = -1.0
C
      If(iireg.ne.0) Return
      econeizol = 0.0
C
C E M C
C
      If(iicel.le.4) Then
         j = iimod
         k = iicel+iitow*4
         Do l=j-1,j+1
            Do m=k-1,k+1
               econeizol = econeizol + emctab(l,m)
            EndDo
         EndDo
         econeizol = 100. * emctab(j,k) / econeizol
C
C H A C 1
C
      ElseIf(iicel.eq.6) Then
         j = iimod
         k = iitow
         Do l=j-1,j+1
            Do m=k-1,k+1
               econeizol = econeizol + hac1tab(l,m)
            EndDo
         EndDo
         econeizol = 100.0 * hac1tab(j,k) / econeizol
C
C H A C 2
C
      ElseIf(iicel.eq.7) Then
         j = iimod
         k = iitow
         Do l=j-1,j+1
            Do m=k-1,k+1
               econeizol = econeizol + hac2tab(l,m)
            EndDo
         EndDo
         econeizol = 100.0 * hac2tab(j,k) / econeizol
C
      EndIf
C
      Return
      End
C


C     ==========================================================
      Real Function TimeWag(Ecell,Izol,Imod,Itow,Icel,Ireg,EImb)
C     ==========================================================
      Implicit None
C
      Integer Imod,Itow,Icel,Ireg
      Real    Ecell,Izol,EImb
C
      Real EMCEcut,HAC1Ecut,HAC2Ecut
      Parameter (EMCEcut=0.5)
      Parameter (HAC1Ecut=1.0)
      Parameter (HAC2Ecut=1.0)
      Real HacIr
      Parameter (HacIr=1.00)
      Real EIzolCut
      Parameter (EIzolCut=51.0)
      Real ImbCut
      Parameter (ImbCut=0.5)
C
      TimeWag = 0.0
C Only FCAL
      If(Ireg.ne.0) Return
      If( (Abs(EImb)/Ecell).gt.ImbCut ) Return
      If(abs(Izol-50.).gt.EIzolCut) Return
      If( ((Icel.le.4).and.(Ecell.gt.EMCEcut)).or.
     +    ((Icel.eq.6).and.(Ecell.gt.HAC1Ecut)).or.
     +    ((Icel.eq.7).and.(Ecell.gt.HAC2Ecut))    ) Then
          TimeWag = 1/(0.40+1.4/(Ecell**0.65))**2
          If((Icel.ge.6).and.(Imod.le.13).and.(Imod.ge.11).and.
     +       (Itow.le.13).and.(Itow.ge.11))
     +       TimeWag = TimeWag * HacIr
      EndIf
C
      Return
      End


C     =======================================
      Real Function getvz(Tf,X,Y,Z,Zint,Tint)
C     =======================================
      Implicit None
C
C     Purpose: for given time (Tf) measured in point (X,Y,Z) and
C              z and time coorditates (Zint,Tint) of nominal vertex
C              it should return coresponded "true" Vertex_z
C
      Real Tf,X,Y,Z,Zint,Tint
C
C Metod description:
C
C    True time of interaction:      Tvtx = (Zint - vz)/c + Tint
C    (vz - seeked vertex)
C
C    Time spend for reaching X,Y,Z: DeltaT = Sqrt(x^2+y^2+(vz-z)^2)/c
C
C    Time shift for cal cells:      ShiftT = - Sqrt(x^2+y^2+z^2)/c
C
C    So we should observe particles in cell in ~
C
C                                   Tf = Tvtx + DeltaT + ShiftT
C    Recalculate we received:
C    ......
C    [c*Tf-c*Tint+Sqrt(x^2+y^2+z^2)-Zint]+vz=Sqrt(x^2+y^2+(vz-z)^2)
C              = A
C     A + vz = Sqrt(x^2+y^2+z^2-2*z*vz+vz^2)   |^2   (A+Vz)>0 ?
C    ......
C     A^2 + 2*A*vz + vz^2 = x^2 + y^2 + z^2 - 2*z*vz + vz^2
C                           \     Ro^2     /
C     2(A+z)*vz = Ro^2 - A^2
C
C     vz = (Ro^2 - A^2)/[2(A+z)]  (this is good if (A+vz)>0 <=> (A+z)>0
C
      Real c
      Parameter(c=30.0)
      Real A,Ro
C
      Ro = Sqrt(x**2+y**2+z**2)
      A  = c*Tf - c*Tint + Ro - Zint
C
      If((A+z).le.0.0) Then
         getvz = 1000000.0
      Else
         getvz = (Ro**2-A**2)/(2*(A+z))
      EndIf
C
      Return
      End



C     ===========================================
      Real Function vtxEfcalcorr(Efcal,FEfcalEMC)
C     ===========================================
      Implicit None
C
      Real Efcal,x,FEfcalEMC,EMCrat
C
      x = log10(Efcal)
      vtxEfcalcorr = -7.9802 + 11.942*x
     +       -5.9376*x**2 + 0.97109*x**3
      vtxEfcalcorr = vtxEfcalcorr +
     +               2.0511 - 4.3378*EMCrat
C
      Return
      End


C     ========================
      Real Function dvz(ncell)
C     ========================
      Implicit None
C
      Integer ncell
      Real x
C
      x = log10(float(ncell))
      If(x.gt.2.) Then
         dvz = 6.3
      Else
         dvz = 29.915 - 24.330*x + 6.2928*x*x
      EndIf
C
      Return
      End
C


C     ============================
      Real Function TfRunCorr(run)
C     ============================
C
      Implicit None
      Integer run,i,inum
      Real rrun,val1,val2,run1,run2
C
      Integer AccRun
      Save AccRun
      Data AccRun / 0 /
      Real    AccShf
      Save AccShf
      Data AccShf / 0.0 /
C
      Integer MaxRCorr94,MaxRCorr95,
     +       MaxRCorr96,MaxRCorr97
      Parameter(MaxRCorr94=     259)
      Parameter(MaxRCorr95=     490)
      Parameter(MaxRCorr96=     679)
      Parameter(MaxRCorr97=     924)
      Integer RunNum94(MaxRCorr94),RunNum95(MaxRCorr95),
     +        RunNum96(MaxRCorr96),RunNum97(MaxRCorr97)
      Real    RunShf94(MaxRCorr94),RunShf95(MaxRCorr95),
     +        RunShf96(MaxRCorr96),RunShf97(MaxRCorr97)
      Integer RunNum94_v1(259),
     +        RunNum95_v1(350),RunNum95_v2(140),
     +        RunNum96_v1(350),RunNum96_v2(329),
     +        RunNum97_v1(350),RunNum97_v2(350),
     +        RunNum97_v3(224)
      Real    RunShf94_v1(259),
     +        RunShf95_v1(350),RunShf95_v2(140),
     +        RunShf96_v1(350),RunShf96_v2(329),
     +        RunShf97_v1(350),RunShf97_v2(350),
     +        RunShf97_v3(224)
      Common /ALALALAM/ RunNum94_v1,
     +                RunNum95_v1,RunNum95_v2,
     +                RunNum96_v1,RunNum96_v2,
     +                RunNum97_v1,RunNum97_v2,RunNum97_v3,
     +                RunShf94_v1,
     +                RunShf95_v1,RunShf95_v2,
     +                RunShf96_v1,RunShf96_v2,
     +                RunShf97_v1,RunShf97_v2,RunShf97_v3
C
      Equivalence (RunShf94,RunShf94_v1)
      Equivalence (RunNum94,RunNum94_v1)
      Equivalence (RunShf95,RunShf95_v1)
      Equivalence (RunNum95,RunNum95_v1)
      Equivalence (RunShf96,RunShf96_v1)
      Equivalence (RunNum96,RunNum96_v1)
      Equivalence (RunShf97,RunShf97_v1)
      Equivalence (RunNum97,RunNum97_v1)
C
C
      Data RunNum94_v1 /
     +      0,   9256,   9272,   9284,   9299,   9301,   9302,
     +   9304,   9306,   9310,   9311,   9316,   9350,   9352,
     +   9360,   9362,   9378,   9380,   9381,   9382,   9385,
     +   9386,   9387,   9393,   9394,   9406,   9407,   9408,
     +   9410,   9413,   9439,   9442,   9456,   9457,   9458,
     +   9466,   9467,   9509,   9511,   9512,   9513,   9519,
     +   9520,   9521,   9526,   9527,   9528,   9536,   9537,
     +   9538,   9543,   9545,   9561,   9566,   9571,   9575,
     +   9576,   9577,   9584,   9586,   9587,   9589,   9591,
     +   9592,   9593,   9606,   9607,   9608,   9609,   9611,
     +   9622,   9623,   9624,   9625,   9626,   9627,   9630,
     +   9631,   9633,   9637,   9638,   9639,   9644,   9645,
     +   9646,   9649,   9651,   9652,   9653,   9654,   9661,
     +   9666,   9667,   9668,   9669,   9679,   9683,   9687,
     +   9688,   9689,   9693,   9694,   9696,   9701,   9702,
     +   9703,   9706,   9708,   9713,   9714,   9715,   9720,
     +   9721,   9722,   9750,   9752,   9753,   9756,   9757,
     +   9758,   9759,   9761,   9767,   9768,   9769,   9773,
     +   9776,   9779,   9783,   9785,   9786,   9792,   9794,
     +   9796,   9800,   9804,   9805,   9806,   9807,   9808,
     +   9812,   9813,   9816,   9817,   9818,   9819,   9823,
     +   9843,   9844,   9845,   9846,   9847,   9848,   9849,
     +   9851,   9852,   9854,   9855,   9856,   9857,   9858,
     +   9861,   9862,   9866,   9867,   9870,   9874,   9890,
     +   9891,   9893,   9894,   9902,   9907,   9908,   9909,
     +   9915,   9916,   9917,   9918,   9920,   9921,   9927,
     +   9960,   9961,   9962,   9965,   9968,   9969,   9976,
     +   9978,   9982,   9988,   9989,   9990,   9991,   9995,
     +   9996,   9997,  10009,  10011,  10012,  10015,  10018,
     +  10019,  10020,  10023,  10029,  10037,  10038,  10049,
     +  10050,  10060,  10062,  10064,  10065,  10073,  10075,
     +  10076,  10078,  10079,  10084,  10085,  10087,  10088,
     +  10092,  10113,  10119,  10134,  10147,  10148,  10149,
     +  10197,  10199,  10201,  10204,  10205,  10208,  10209,
     +  10210,  10211,  10212,  10215,  10220,  10229,  10231,
     +  10235,  10236,  10242,  10243,  10244,  10245,  10246,
     +  10249,  10259,  10260,  10261,  10262,  10263,  10328 /
      Data RunNum95_v1 /
     +  10330,  11541,  11543,  11547,  11548,  11550,  11562,
     +  11563,  11571,  11573,  11586,  11587,  11592,  11593,
     +  11594,  11595,  11600,  11635,  11636,  11640,  11642,
     +  11645,  11646,  11661,  11663,  11670,  11675,  11682,
     +  11683,  11684,  11711,  11713,  11714,  11715,  11716,
     +  11805,  11806,  11807,  11815,  11816,  11817,  11818,
     +  11819,  11820,  11821,  11822,  11824,  11825,  11826,
     +  11841,  12012,  12018,  12019,  12020,  12021,  12025,
     +  12026,  12049,  12075,  12076,  12077,  12119,  12120,
     +  12131,  12139,  12140,  12143,  12144,  12146,  12148,
     +  12157,  12158,  12159,  12161,  12162,  12163,  12164,
     +  12171,  12172,  12173,  12187,  12188,  12190,  12192,
     +  12199,  12208,  12209,  12211,  12213,  12214,  12215,
     +  12232,  12244,  12246,  12247,  12248,  12249,  12256,
     +  12262,  12265,  12267,  12268,  12300,  12302,  12313,
     +  12314,  12331,  12332,  12334,  12365,  12366,  12370,
     +  12418,  12424,  12425,  12426,  12430,  12432,  12447,
     +  12451,  12452,  12453,  12454,  12456,  12469,  12470,
     +  12471,  12486,  12487,  12488,  12489,  12493,  12495,
     +  12496,  12497,  12508,  12509,  12510,  12513,  12514,
     +  12537,  12538,  12539,  12541,  12546,  12547,  12548,
     +  12561,  12562,  12563,  12564,  12571,  12576,  12613,
     +  12614,  12616,  12619,  12622,  12629,  12630,  12631,
     +  12635,  12636,  12637,  12640,  12641,  12642,  12650,
     +  12671,  12672,  12673,  12682,  12683,  12685,  12693,
     +  12695,  12696,  12697,  12704,  12705,  12706,  12707,
     +  12717,  12718,  12729,  12732,  12733,  12734,  12737,
     +  12738,  12748,  12749,  12750,  12751,  12753,  12754,
     +  12758,  12759,  12760,  12761,  12788,  12789,  12793,
     +  12794,  12795,  12808,  12825,  12827,  12828,  12829,
     +  12830,  12831,  12853,  12854,  12855,  12872,  12874,
     +  12884,  12885,  12886,  12895,  12896,  12897,  12905,
     +  12906,  12907,  12908,  12909,  12910,  12956,  12957,
     +  12969,  12973,  12982,  12983,  12985,  12990,  12991,
     +  12998,  12999,  13002,  13004,  13007,  13008,  13020,
     +  13035,  13045,  13046,  13047,  13048,  13051,  13057,
     +  13085,  13086,  13088,  13093,  13100,  13101,  13104,
     +  13108,  13109,  13114,  13115,  13116,  13118,  13119,
     +  13121,  13122,  13123,  13126,  13130,  13134,  13135,
     +  13136,  13144,  13145,  13146,  13149,  13150,  13151,
     +  13152,  13153,  13156,  13157,  13158,  13171,  13173,
     +  13175,  13176,  13179,  13180,  13182,  13183,  13186,
     +  13187,  13200,  13201,  13202,  13203,  13213,  13214,
     +  13215,  13216,  13217,  13218,  13223,  13224,  13229,
     +  13230,  13237,  13239,  13241,  13242,  13248,  13254,
     +  13255,  13256,  13271,  13284,  13285,  13288,  13297,
     +  13298,  13299,  13308,  13315,  13316,  13319,  13327,
     +  13328,  13329,  13335,  13336,  13348,  13349,  13350,
     +  13351,  13370,  13371,  13372,  13375,  13376,  13390,
     +  13391,  13392,  13393,  13397,  13399,  13400,  13401 /
      Data RunNum95_v2 /
     +  13402,  13404,  13409,  13410,  13420,  13421,  13424,
     +  13437,  13438,  13440,  13441,  13445,  13458,  13459,
     +  13461,  13463,  13472,  13475,  13478,  13479,  13480,
     +  13487,  13488,  13492,  13496,  13501,  13527,  13551,
     +  13556,  13558,  13564,  13573,  13582,  13598,  13599,
     +  13601,  13604,  13605,  13612,  13654,  13655,  13657,
     +  13658,  13659,  13661,  13663,  13673,  13674,  13675,
     +  13677,  13678,  13690,  13699,  13702,  13704,  13726,
     +  13732,  13735,  13736,  13741,  13749,  13750,  13751,
     +  13752,  13761,  13762,  13768,  13769,  13775,  13778,
     +  13779,  13794,  13795,  13796,  13802,  13804,  13805,
     +  13807,  13817,  13818,  13839,  13844,  13852,  13854,
     +  13875,  13876,  13878,  13879,  13880,  13881,  13889,
     +  13890,  13908,  13909,  13910,  13911,  13912,  13913,
     +  13915,  13950,  13952,  13959,  13968,  13984,  13990,
     +  13991,  13993,  14018,  14027,  14028,  14035,  14036,
     +  14048,  14054,  14056,  14219,  14220,  14234,  14235,
     +  14236,  14244,  14245,  14246,  14248,  14253,  14254,
     +  14264,  14265,  14270,  14304,  14361,  14362,  14375,
     +  14376,  14394,  14395,  19996,  19997,  19998,  19999 /
      Data RunNum96_v1 /
     +  20000,  20718,  20719,  20721,  20727,  20728,  20735,
     +  20750,  20751,  20752,  20753,  20754,  20755,  20771,
     +  20785,  20786,  20787,  20790,  20791,  20792,  20819,
     +  20820,  20857,  20858,  20864,  20869,  20870,  20871,
     +  20872,  20875,  20876,  20879,  20905,  20907,  20910,
     +  20911,  21165,  21166,  21167,  21168,  21186,  21187,
     +  21189,  21190,  21191,  21192,  21193,  21202,  21203,
     +  21204,  21205,  21206,  21207,  21209,  21217,  21218,
     +  21219,  21226,  21227,  21260,  21261,  21262,  21263,
     +  21278,  21279,  21285,  21286,  21288,  21289,  21290,
     +  21296,  21297,  21298,  21303,  21304,  21338,  21339,
     +  21356,  21357,  21371,  21373,  21374,  21375,  21376,
     +  21377,  21378,  21383,  21385,  21386,  21392,  21393,
     +  21440,  21443,  21445,  21447,  21448,  21449,  21469,
     +  21470,  21471,  21472,  21473,  21477,  21478,  21495,
     +  21496,  21497,  21504,  21505,  21506,  21508,  21512,
     +  21513,  21514,  21516,  21521,  21523,  21524,  21525,
     +  21526,  21528,  21529,  21545,  21546,  21548,  21551,
     +  21553,  21555,  21556,  21558,  21559,  21560,  21561,
     +  21562,  21566,  21567,  21569,  21570,  21576,  21577,
     +  21578,  21579,  21580,  21581,  21582,  21587,  21588,
     +  21589,  21590,  21591,  21609,  21610,  21612,  21613,
     +  21614,  21615,  21616,  21617,  21618,  21626,  21629,
     +  21630,  21631,  21634,  21635,  21636,  21637,  21639,
     +  21640,  21641,  21642,  21646,  21647,  21648,  21649,
     +  21652,  21653,  21654,  21659,  21660,  21662,  21663,
     +  21665,  21669,  21670,  21689,  21690,  21702,  21703,
     +  21704,  21705,  21706,  21707,  21708,  21709,  21713,
     +  21714,  21716,  21717,  21718,  21719,  21720,  21724,
     +  21725,  21726,  21727,  21728,  21732,  21733,  21734,
     +  21735,  21736,  21737,  21745,  21746,  21752,  21755,
     +  21756,  21757,  21758,  21759,  21763,  21767,  21777,
     +  21778,  21784,  21785,  21786,  21791,  21792,  21793,
     +  21794,  21795,  21796,  21797,  21798,  21813,  21817,
     +  21819,  21821,  21838,  21839,  21840,  21841,  21848,
     +  21849,  21850,  21851,  21853,  21871,  21872,  21873,
     +  21874,  21876,  21878,  21883,  21884,  21885,  21886,
     +  21887,  21890,  21892,  21893,  21894,  21895,  21896,
     +  21897,  21898,  21899,  21900,  21901,  21904,  21905,
     +  21906,  21907,  21910,  21911,  21912,  21913,  21914,
     +  21915,  21928,  21929,  21930,  21933,  21938,  21939,
     +  21940,  21941,  21942,  21943,  21944,  21945,  21949,
     +  21950,  21951,  21953,  21957,  21958,  21960,  21964,
     +  21967,  21973,  21974,  21975,  21976,  21979,  21989,
     +  21990,  22000,  22001,  22002,  22003,  22005,  22018,
     +  22019,  22020,  22025,  22027,  22028,  22029,  22030,
     +  22031,  22032,  22036,  22045,  22047,  22048,  22050,
     +  22051,  22052,  22058,  22059,  22061,  22063,  22065,
     +  22066,  22067,  22072,  22075,  22076,  22077,  22078,
     +  22079,  22080,  22087,  22089,  22090,  22091,  22092 /
      Data RunNum96_v2 /
     +  22094,  22111,  22112,  22113,  22114,  22115,  22116,
     +  22124,  22125,  22126,  22127,  22128,  22129,  22130,
     +  22136,  22137,  22140,  22141,  22142,  22147,  22148,
     +  22149,  22150,  22151,  22152,  22153,  22156,  22169,
     +  22170,  22171,  22172,  22173,  22174,  22186,  22187,
     +  22188,  22189,  22192,  22204,  22206,  22207,  22208,
     +  22209,  22211,  22215,  22216,  22217,  22218,  22219,
     +  22220,  22221,  22224,  22225,  22226,  22227,  22233,
     +  22234,  22235,  22237,  22238,  22239,  22240,  22241,
     +  22242,  22244,  22245,  22246,  22247,  22248,  22250,
     +  22251,  22252,  22253,  22254,  22260,  22261,  22262,
     +  22263,  22264,  22265,  22266,  22269,  22270,  22273,
     +  22274,  22275,  22276,  22277,  22278,  22279,  22280,
     +  22281,  22289,  22290,  22291,  22292,  22294,  22295,
     +  22297,  22300,  22301,  22302,  22310,  22311,  22312,
     +  22336,  22337,  22338,  22339,  22340,  22345,  22346,
     +  22347,  22353,  22354,  22362,  22364,  22366,  22367,
     +  22368,  22403,  22407,  22408,  22409,  22414,  22415,
     +  22416,  22417,  22419,  22420,  22423,  22424,  22425,
     +  22426,  22437,  22438,  22439,  22441,  22442,  22447,
     +  22451,  22452,  22454,  22458,  22459,  22460,  22462,
     +  22466,  22468,  22469,  22479,  22490,  22491,  22492,
     +  22493,  22494,  22495,  22496,  22505,  22513,  22514,
     +  22538,  22539,  22548,  22549,  22555,  22556,  22557,
     +  22559,  22562,  22563,  22571,  22574,  22575,  22576,
     +  22577,  22583,  22584,  22596,  22597,  22598,  22599,
     +  22600,  22601,  22605,  22606,  22610,  22611,  22612,
     +  22613,  22614,  22615,  22617,  22624,  22625,  22627,
     +  22636,  22639,  22640,  22647,  22648,  22649,  22650,
     +  22651,  22652,  22653,  22654,  22655,  22656,  22657,
     +  22658,  22660,  22675,  22679,  22680,  22681,  22682,
     +  22685,  22686,  22695,  22696,  22698,  22704,  22708,
     +  22709,  22710,  22711,  22713,  22715,  22716,  22717,
     +  22718,  22719,  22722,  22723,  22724,  22725,  22726,
     +  22731,  22732,  22736,  22737,  22740,  22745,  22746,
     +  22747,  22753,  22756,  22757,  22764,  22765,  22767,
     +  22768,  22769,  22777,  22780,  22795,  22796,  22810,
     +  22812,  22813,  22814,  22818,  22819,  22820,  22821,
     +  22822,  22824,  22835,  22836,  22837,  22838,  22839,
     +  22840,  22841,  22843,  22848,  22849,  22850,  22851,
     +  22856,  22858,  22860,  22861,  22862,  22863,  22864,
     +  22868,  22869,  22871,  22872,  22873,  22874,  22879,
     +  22880,  22881,  22882,  22883,  22884,  22885,  22886,
     +  22888,  22890,  22892,  22894,  22898,  22901,  22902,
     +  22903,  22904,  22906,  22908,  22909,  22916,  22917,
     +  22918,  22920,  22921,  22922,  22934,  22942,  22943,
     +  22944,  22945,  22946,  22947,  22953,  22954,  24999 /
      Data RunNum97_v1 /
     +  25000,  25190,  25196,  25197,  25198,  25200,  25202,
     +  25204,  25205,  25211,  25212,  25214,  25220,  25224,
     +  25225,  25237,  25238,  25243,  25244,  25246,  25247,
     +  25269,  25270,  25271,  25279,  25280,  25281,  25282,
     +  25291,  25292,  25304,  25312,  25316,  25317,  25319,
     +  25320,  25346,  25352,  25353,  25354,  25360,  25368,
     +  25371,  25372,  25377,  25378,  25379,  25380,  25386,
     +  25387,  25388,  25393,  25394,  25396,  25397,  25400,
     +  25403,  25404,  25405,  25410,  25416,  25419,  25420,
     +  25421,  25431,  25432,  25437,  25438,  25440,  25441,
     +  25443,  25448,  25456,  25459,  25461,  25464,  25472,
     +  25473,  25474,  25475,  25480,  25481,  25482,  25483,
     +  25484,  25487,  25492,  25493,  25494,  25497,  25515,
     +  25516,  25518,  25525,  25526,  25541,  25542,  25544,
     +  25555,  25561,  25562,  25563,  25576,  25577,  25578,
     +  25579,  25580,  25596,  25597,  25598,  25612,  25613,
     +  25625,  25626,  25634,  25640,  25641,  25642,  25645,
     +  25646,  25653,  25655,  25656,  25657,  25658,  25659,
     +  25660,  25661,  25662,  25686,  25693,  25694,  25695,
     +  25696,  25697,  25698,  25699,  25713,  25718,  25719,
     +  25721,  25722,  25723,  25724,  25725,  25730,  25731,
     +  25732,  25741,  25745,  25746,  25747,  25748,  25749,
     +  25753,  25760,  25761,  25762,  25764,  25765,  25766,
     +  25767,  25784,  25785,  25786,  25787,  25788,  25789,
     +  25790,  25792,  25804,  25805,  25808,  25809,  25810,
     +  25821,  25831,  25832,  25833,  25834,  25835,  25839,
     +  25840,  25841,  25842,  25844,  25846,  25847,  25852,
     +  25854,  25855,  25861,  25862,  25863,  25864,  25866,
     +  25867,  25872,  25873,  25874,  25877,  25878,  25879,
     +  25882,  25884,  25885,  25886,  25887,  25894,  25895,
     +  25896,  25897,  25898,  25899,  25903,  25906,  25916,
     +  25917,  25918,  25923,  25925,  25926,  25934,  25935,
     +  25936,  25939,  25940,  25941,  25946,  25947,  25948,
     +  25949,  25952,  25955,  25956,  25957,  25958,  25963,
     +  25968,  25972,  25973,  25977,  25978,  25979,  25980,
     +  25981,  25992,  25993,  25994,  25997,  26002,  26011,
     +  26017,  26018,  26027,  26034,  26036,  26043,  26046,
     +  26047,  26052,  26053,  26054,  26055,  26060,  26061,
     +  26062,  26063,  26065,  26066,  26067,  26071,  26072,
     +  26073,  26074,  26086,  26097,  26098,  26099,  26100,
     +  26102,  26109,  26112,  26113,  26114,  26115,  26117,
     +  26118,  26122,  26123,  26125,  26128,  26129,  26131,
     +  26132,  26133,  26137,  26138,  26139,  26145,  26146,
     +  26147,  26148,  26150,  26155,  26168,  26169,  26177,
     +  26178,  26179,  26180,  26181,  26185,  26186,  26198,
     +  26206,  26207,  26208,  26215,  26217,  26218,  26221,
     +  26227,  26228,  26229,  26230,  26231,  26233,  26234,
     +  26244,  26245,  26246,  26271,  26272,  26273,  26275,
     +  26276,  26277,  26281,  26282,  26284,  26285,  26289,
     +  26290,  26291,  26292,  26299,  26300,  26302,  26307 /
      Data RunNum97_v2 /
     +  26308,  26310,  26311,  26312,  26313,  26314,  26320,
     +  26322,  26324,  26325,  26326,  26330,  26331,  26332,
     +  26333,  26334,  26343,  26345,  26346,  26347,  26357,
     +  26358,  26359,  26360,  26362,  26363,  26364,  26365,
     +  26370,  26371,  26372,  26380,  26383,  26384,  26387,
     +  26393,  26396,  26397,  26398,  26399,  26400,  26401,
     +  26402,  26403,  26411,  26412,  26413,  26414,  26415,
     +  26420,  26422,  26432,  26434,  26435,  26437,  26440,
     +  26452,  26453,  26454,  26455,  26456,  26457,  26458,
     +  26460,  26461,  26462,  26463,  26464,  26471,  26472,
     +  26473,  26474,  26475,  26476,  26486,  26487,  26488,
     +  26489,  26490,  26491,  26492,  26493,  26494,  26495,
     +  26497,  26500,  26501,  26504,  26505,  26506,  26513,
     +  26514,  26515,  26524,  26525,  26529,  26531,  26536,
     +  26537,  26538,  26539,  26543,  26544,  26545,  26549,
     +  26550,  26551,  26552,  26553,  26554,  26558,  26560,
     +  26568,  26570,  26571,  26572,  26666,  26667,  26671,
     +  26672,  26673,  26674,  26675,  26680,  26681,  26682,
     +  26683,  26685,  26686,  26687,  26688,  26690,  26691,
     +  26701,  26702,  26705,  26707,  26709,  26711,  26712,
     +  26714,  26720,  26722,  26728,  26733,  26734,  26737,
     +  26744,  26745,  26752,  26757,  26762,  26763,  26764,
     +  26765,  26769,  26774,  26775,  26776,  26784,  26787,
     +  26789,  26793,  26794,  26795,  26796,  26800,  26801,
     +  26818,  26819,  26823,  26828,  26829,  26830,  26831,
     +  26832,  26837,  26838,  26840,  26841,  26842,  26843,
     +  26864,  26865,  26871,  26872,  26873,  26882,  26884,
     +  26885,  26886,  26890,  26892,  26893,  26900,  26904,
     +  26905,  26907,  26908,  26909,  26910,  26911,  26912,
     +  26924,  26925,  26933,  26940,  26941,  26942,  26945,
     +  26946,  26947,  26949,  26959,  26960,  26966,  26967,
     +  26968,  26969,  26970,  26971,  26976,  26978,  26979,
     +  26980,  26981,  26983,  26987,  26990,  26991,  26992,
     +  26993,  27003,  27017,  27018,  27019,  27020,  27021,
     +  27026,  27028,  27029,  27030,  27034,  27035,  27047,
     +  27048,  27049,  27050,  27051,  27052,  27053,  27066,
     +  27067,  27068,  27076,  27077,  27078,  27081,  27090,
     +  27091,  27092,  27093,  27094,  27096,  27097,  27100,
     +  27101,  27105,  27134,  27135,  27137,  27138,  27142,
     +  27143,  27147,  27148,  27149,  27152,  27162,  27163,
     +  27164,  27166,  27167,  27173,  27177,  27178,  27179,
     +  27185,  27186,  27187,  27188,  27190,  27191,  27228,
     +  27229,  27230,  27231,  27232,  27233,  27234,  27236,
     +  27240,  27241,  27242,  27243,  27265,  27266,  27267,
     +  27268,  27269,  27272,  27273,  27274,  27291,  27305,
     +  27306,  27307,  27312,  27313,  27314,  27318,  27320,
     +  27321,  27322,  27327,  27336,  27337,  27338,  27339,
     +  27340,  27345,  27346,  27347,  27348,  27349,  27350,
     +  27351,  27352,  27353,  27354,  27355,  27356,  27362,
     +  27363,  27364,  27367,  27373,  27374,  27377,  27378 /
      Data RunNum97_v3 /
     +  27379,  27383,  27385,  27386,  27387,  27388,  27400,
     +  27401,  27402,  27407,  27408,  27409,  27412,  27413,
     +  27415,  27419,  27420,  27421,  27425,  27426,  27427,
     +  27428,  27429,  27430,  27435,  27436,  27440,  27441,
     +  27445,  27446,  27450,  27451,  27454,  27455,  27456,
     +  27461,  27462,  27463,  27470,  27471,  27475,  27476,
     +  27477,  27481,  27482,  27485,  27488,  27489,  27490,
     +  27507,  27509,  27525,  27526,  27527,  27528,  27529,
     +  27534,  27535,  27538,  27539,  27540,  27553,  27554,
     +  27555,  27559,  27560,  27561,  27562,  27563,  27564,
     +  27565,  27570,  27571,  27572,  27573,  27574,  27577,
     +  27578,  27580,  27581,  27585,  27589,  27590,  27592,
     +  27593,  27594,  27596,  27601,  27603,  27610,  27618,
     +  27619,  27620,  27621,  27622,  27623,  27624,  27625,
     +  27626,  27630,  27631,  27632,  27635,  27640,  27641,
     +  27642,  27643,  27644,  27645,  27646,  27647,  27651,
     +  27652,  27653,  27654,  27655,  27656,  27660,  27661,
     +  27662,  27663,  27668,  27669,  27670,  27676,  27677,
     +  27684,  27685,  27686,  27687,  27697,  27698,  27703,
     +  27704,  27705,  27706,  27708,  27713,  27716,  27717,
     +  27718,  27719,  27720,  27723,  27725,  27726,  27727,
     +  27728,  27729,  27730,  27731,  27732,  27737,  27738,
     +  27739,  27740,  27741,  27743,  27744,  27745,  27746,
     +  27761,  27762,  27763,  27766,  27767,  27769,  27770,
     +  27773,  27774,  27781,  27782,  27783,  27785,  27786,
     +  27793,  27794,  27795,  27796,  27797,  27798,  27799,
     +  27804,  27805,  27812,  27813,  27814,  27815,  27816,
     +  27818,  27819,  27823,  27824,  27826,  27827,  27828,
     +  27829,  27830,  27831,  27835,  27836,  27837,  27838,
     +  27839,  27862,  27863,  27866,  27867,  27868,  27870,
     +  27871,  27875,  27876,  27877,  27878,  27879,  27880,
     +  27881,  28995,  28996,  28997,  28998,  28999,  29000 /
      Data RunShf94_v1 /
     +  0.000, -1.023, -1.023, -0.855, -0.855, -0.701, -0.701,
     + -1.307, -1.307, -1.322, -1.038, -1.038, -0.268, -0.268,
     + -0.944, -0.944, -0.413, -0.413, -0.847, -0.399, -0.399,
     + -0.583, -0.797, -0.784, -0.784, -1.056, -1.056, -1.072,
     + -1.258, -1.258, -1.113, -1.175, -1.175, -0.872, -1.155,
     + -1.155, -1.420, -0.391, -0.391, -0.471, -0.372, -0.372,
     + -0.537, -0.537, -0.861, -0.861, -0.373, -0.373, -0.674,
     + -0.394, -0.394, -0.443, -0.443, -0.429, -0.429, -0.728,
     + -0.728, -0.453, -0.453, -0.348, -0.702, -0.667, -0.667,
     + -0.446, -0.269, -0.269, -0.201, -0.201, -0.384, -0.384,
     + -0.396, -0.396, -0.712, -0.231, -0.231, -0.556, -0.556,
     + -0.297, -0.548, -0.548, -0.555, -0.555, -0.111, -0.317,
     + -0.360, -0.408, -0.408,  0.110,  0.110, -0.584, -0.584,
     + -0.347, -0.347, -0.495, -0.495, -0.519, -0.519, -0.519,
     + -0.379, -0.190, -0.190,  0.056,  0.056, -0.348, -0.348,
     + -0.327, -0.327, -0.556, -0.556, -0.578, -0.756, -0.756,
     + -0.650, -0.650, -0.442, -0.577, -0.577, -1.023, -0.719,
     + -0.465, -0.682, -0.614, -0.604, -0.378, -0.950, -0.733,
     + -0.733, -0.744, -0.547, -0.774, -0.930, -0.930, -0.704,
     + -0.704, -0.589, -0.529, -0.795, -0.713, -0.713, -0.654,
     + -0.528, -0.528, -0.219, -0.219, -0.480, -0.595, -0.595,
     + -0.643, -0.643, -0.733, -0.733, -0.339, -0.339, -0.576,
     + -0.606, -0.536, -0.540, -0.540, -0.574, -0.948, -0.444,
     + -0.444, -0.399, -0.506, -0.745, -0.614, -0.614, -0.564,
     + -0.684, -0.684, -0.351, -0.622, -0.622, -0.597, -0.775,
     + -0.775, -0.874, -0.874, -0.438, -0.438, -0.832, -0.504,
     + -0.650, -0.650, -0.515, -0.575, -0.559, -0.594, -0.625,
     + -0.625, -0.528, -0.756, -0.699, -0.694, -0.700, -0.700,
     + -0.533, -0.533, -0.859, -0.881, -0.881, -0.991, -0.959,
     + -0.763, -0.529, -0.806, -0.806, -0.672, -1.034, -1.034,
     + -0.603, -0.661, -0.792, -0.792, -0.810, -0.810, -0.833,
     + -0.649, -0.649, -0.842, -0.842, -0.800, -0.800, -0.852,
     + -0.852, -1.088, -1.088, -0.695, -0.819, -0.388, -0.891,
     + -0.891, -0.809, -0.933, -0.933, -0.864, -0.878, -0.878,
     + -1.131, -1.131, -1.664, -1.041, -0.940, -0.940, -1.021,
     + -1.021, -1.104, -0.739, -0.868, -0.709, -1.035, -0.888,
     + -0.931, -0.931, -1.075, -1.075, -1.185, -1.185,  0.000 /
      Data RunShf95_v1 /
     +  0.000,  0.828,  0.828,  0.591,  0.502,  1.110,  1.110,
     +  0.756,  0.756,  0.743,  0.743,  0.501,  0.501,  0.580,
     +  0.658,  0.658,  0.480,  0.480,  0.224,  0.224,  0.594,
     +  0.594,  0.705,  0.705,  0.249,  0.249,  0.593,  0.593,
     +  0.974,  0.974,  0.291,  0.291,  0.634,  0.569,  0.494,
     +  0.494,  0.622,  0.531,  0.531,  0.824,  0.824,  0.723,
     +  0.452,  0.452,  0.392,  0.392,  0.627,  0.627,  0.161,
     +  0.161,  0.511,  0.511,  0.324,  0.324,  0.549,  0.549,
     +  0.429,  0.429,  0.628,  0.570,  0.383,  0.383,  0.746,
     +  0.746,  0.310,  0.404,  0.404,  0.206,  0.523,  0.523,
     +  0.439,  0.439,  0.567,  0.345,  0.345,  0.582,  0.359,
     +  0.359,  0.405,  0.434,  0.432,  0.672,  0.949,  0.359,
     +  0.359,  0.730,  0.730,  0.493,  0.493,  0.788,  0.788,
     +  0.575,  0.575,  0.287,  0.580,  0.580,  0.553,  0.553,
     +  0.526,  0.648,  0.648,  0.526,  0.430,  0.084,  0.695,
     +  0.624,  0.624,  0.313,  0.725,  0.725,  0.519,  0.519,
     +  0.379,  0.379,  0.425,  0.546,  0.546,  0.499,  0.483,
     +  0.483,  0.542,  0.542,  0.517,  0.584,  0.584,  0.408,
     +  0.372,  0.448,  0.448,  0.166,  0.166,  0.695,  0.329,
     +  0.275,  0.275,  0.236,  0.342,  0.258, -0.094, -0.094,
     +  0.369,  0.369,  0.552,  0.552,  0.402,  0.582,  0.335,
     +  0.467,  0.423,  0.258,  0.460,  0.460,  0.630,  0.630,
     +  0.611,  0.579,  0.371,  0.391,  0.362,  0.149,  0.539,
     +  0.539,  0.257,  0.278,  0.421,  0.530,  0.297,  0.297,
     +  0.256,  0.294,  0.294,  0.324,  0.353,  0.353,  0.378,
     +  0.279,  0.279,  0.166,  0.278,  0.262,  0.361,  0.361,
     +  0.506,  0.506,  0.573,  0.310,  0.470,  0.249,  0.527,
     +  0.280,  0.550,  0.475,  0.672,  0.619,  0.619,  0.274,
     +  0.376,  0.376,  0.653,  0.496,  0.678,  0.678,  0.681,
     +  0.532,  0.444,  0.444,  0.510,  0.510,  0.628,  0.515,
     +  0.602,  0.602,  0.496,  0.555,  0.807,  0.576,  0.477,
     +  0.477,  0.430,  0.664,  0.767,  1.352,  1.352,  0.704,
     +  0.677,  0.329,  0.689,  0.618,  0.563,  0.563,  0.570,
     +  0.436,  0.558,  0.558,  0.283,  0.418,  0.418,  0.701,
     +  0.701,  0.381,  0.381,  0.689,  0.414,  0.414,  0.471,
     +  0.696,  0.696,  0.362,  0.362,  0.563,  0.563,  0.485,
     +  0.456,  0.456,  0.836,  0.836,  0.283,  0.727,  0.685,
     +  0.685,  0.794,  0.597,  0.620,  0.620,  0.444,  0.416,
     +  0.796,  0.796,  0.824,  0.397,  0.469,  0.469,  0.342,
     +  0.460,  0.554,  0.725,  0.725,  0.583,  0.583,  0.440,
     +  0.563,  0.563,  0.716,  0.716,  0.666,  0.666,  0.694,
     +  0.694,  0.489,  0.489,  0.503,  0.503,  0.388,  0.388,
     +  0.551,  0.551,  0.238,  0.238,  0.670,  0.666,  0.834,
     +  0.834,  0.463,  0.463,  0.447,  0.692,  0.731,  0.132,
     +  0.451,  0.451,  0.337,  0.513,  0.513,  0.419,  0.419,
     +  0.587,  0.567,  0.567,  0.392,  0.392,  0.220,  0.220,
     +  0.432,  0.674,  0.628,  0.530,  0.325,  0.117,  0.774,
     +  0.620,  0.620,  0.466,  0.883,  0.883,  0.709,  0.925,
     +  0.254,  0.254,  0.749,  0.749,  0.791,  0.066,  0.066,
     +  0.076,  0.076,  0.569,  0.569,  0.736,  0.736,  0.816 /
      Data RunShf95_v2 /
     +  0.816,  0.201,  0.201,  0.335,  0.660,  0.660,  0.529,
     +  0.432,  0.432,  0.650,  0.527,  1.097,  0.577,  0.709,
     +  0.625,  0.625,  0.637,  0.524,  0.328,  0.509,  0.729,
     +  0.729,  0.616,  0.642,  0.786,  0.686,  0.644,  0.732,
     +  0.732,  0.840,  0.840,  0.544,  0.544,  0.659,  0.659,
     +  0.398,  0.999,  0.999,  0.488,  0.390,  0.390,  0.733,
     +  0.753,  0.753,  0.905,  0.905,  2.011,  2.011,  2.039,
     +  2.039,  0.877,  0.877,  0.611,  0.312,  0.502,  0.716,
     +  0.765,  0.543,  0.635,  0.640,  0.461,  0.461,  0.500,
     +  0.477,  0.477,  0.273,  0.339,  0.477,  0.567,  0.567,
     +  0.540,  0.476,  0.476,  0.294,  0.327,  0.327,  0.496,
     +  0.459,  0.650,  0.710,  0.710,  0.343,  0.380,  0.600,
     +  0.458,  0.299,  0.299,  0.447,  0.447,  0.319,  0.247,
     +  0.247,  0.088,  0.593,  0.235,  0.501,  0.209,  0.315,
     +  0.397,  0.397,  0.184,  0.351,  0.351,  0.518,  0.518,
     +  0.172,  0.458,  0.316,  0.316,  0.610,  0.294,  0.245,
     +  0.519,  0.514,  0.449,  0.449,  0.609,  0.609,  0.337,
     +  0.775,  0.501,  0.501,  0.445,  0.445,  0.366,  0.366,
     +  1.327,  1.327,  0.231,  0.231,  0.440,  0.717,  0.717,
     +  0.537,  0.537,  0.643,  0.000,  0.000,  0.000,  0.000 /
      Data RunShf96_v1 /
     +  0.000,  0.168,  0.168,  0.239,  0.239,  0.325,  0.325,
     +  0.112,  0.486,  0.486,  0.133,  0.416,  0.390,  0.125,
     +  0.125,  0.020,  0.618,  0.618, -0.039,  0.032,  0.032,
     +  0.018,  0.018, -0.055, -0.055, -0.049, -0.409, -0.409,
     +  0.085,  0.085,  0.065,  0.065, -0.072, -0.225, -0.225,
     + -0.465, -0.465, -0.298, -0.310, -0.310, -0.091, -0.091,
     + -0.339, -0.325, -0.325, -0.318, -0.529, -0.222, -0.222,
     + -0.286, -0.421, -0.421, -0.574, -0.574, -0.228, -0.017,
     + -0.017, -0.193, -0.193,  0.307,  0.307,  0.255,  0.255,
     +  0.489,  0.489, -0.161,  0.059,  0.059,  0.373,  0.624,
     +  1.126,  1.126,  0.955,  0.955,  0.995,  0.995,  0.342,
     +  0.209,  0.344,  0.344,  0.165,  0.165,  0.376,  0.729,
     +  0.223,  0.399,  0.361,  0.361,  0.124,  0.172,  0.162,
     + -0.367, -0.367, -0.134, -0.134,  0.190,  0.192,  0.192,
     + -0.282, -0.282,  0.113,  0.113,  0.064,  0.064,  0.144,
     +  0.144,  0.043,  0.043,  0.101,  0.153,  0.153,  0.415,
     + -0.037, -0.037, -0.019, -0.213, -0.213, -0.237,  0.085,
     + -0.232, -0.232, -0.309, -0.309,  0.011,  0.011, -0.341,
     + -0.341, -0.091, -0.091, -0.118, -0.118, -0.875, -0.875,
     + -0.129, -0.129, -0.106,  0.083,  0.083,  0.209, -0.102,
     + -0.102, -0.075, -0.075, -0.207, -0.207,  0.023,  0.023,
     + -0.195, -0.195,  0.158, -0.181, -0.013, -0.013, -0.339,
     +  0.061,  0.049,  0.049, -0.316, -0.316, -0.095, -0.095,
     + -0.008, -0.008, -0.289, -0.331, -0.331,  0.029, -0.129,
     + -0.150, -0.150, -0.144, -0.271, -0.178, -0.045, -0.338,
     + -0.198, -0.122, -0.130, -0.130, -0.068, -0.202, -0.220,
     + -0.205, -0.145, -0.063, -0.063, -0.115,  0.021, -0.146,
     + -0.219, -0.219, -0.183, -0.128, -0.128, -0.151, -0.151,
     + -0.093, -0.424, -0.135, -0.263, -0.263, -0.333, -0.184,
     + -0.333, -0.270, -0.270, -0.121,  0.037,  0.037, -0.031,
     +  0.073,  0.075,  0.107,  0.107,  0.106,  0.023,  0.023,
     + -0.006, -0.006,  0.204,  0.052,  0.321,  0.016,  0.016,
     +  0.522,  0.522,  0.448,  0.448,  0.592,  0.313,  0.313,
     +  0.256,  0.298,  0.431,  0.139,  0.139,  0.282,  0.282,
     +  0.315,  0.315,  0.320,  0.320,  0.306,  0.433,  0.433,
     +  0.433,  0.222,  0.222,  0.458,  0.442,  0.290,  0.130,
     +  0.383,  0.383,  0.302,  0.343,  0.267,  0.039,  0.216,
     +  0.326,  0.326,  0.458,  0.458, -0.016,  0.358,  0.358,
     +  0.294,  0.294,  0.259,  0.259,  0.230,  0.225,  0.163,
     +  0.185,  0.682,  0.682, -0.112,  0.177,  0.273,  0.205,
     +  0.205,  0.149,  0.149,  0.301,  0.291,  0.140,  0.140,
     +  0.239,  0.343,  0.343,  0.360,  0.360,  0.169,  0.169,
     +  0.419,  0.419,  0.221,  0.586,  0.650,  0.664,  0.664,
     +  0.783,  0.529,  0.597,  0.952,  0.952,  0.725,  0.548,
     +  0.553,  0.688,  0.688,  0.523,  0.835,  1.006,  0.707,
     +  0.589,  0.676,  0.676,  0.553,  0.472,  0.573,  0.573,
     +  0.769,  1.045,  0.782,  0.782,  0.822,  0.822,  0.980,
     +  0.697,  0.697,  0.588,  0.715,  0.569,  0.355,  0.668,
     +  0.579,  0.579,  0.767,  0.767,  0.565,  0.565,  0.264,
     +  0.630,  0.632,  0.632,  0.722,  0.722,  0.759,  0.759 /
      Data RunShf96_v2 /
     +  0.674,  0.608,  0.675,  0.435,  0.673,  0.613,  0.613,
     +  0.635,  0.635,  0.773,  0.773,  0.806,  0.713,  0.713,
     +  0.423,  0.530,  0.530,  0.505,  0.283,  0.283,  0.629,
     +  0.421,  0.421,  0.577,  0.685,  0.753,  0.753, -1.019,
     + -0.876, -0.379, -0.398, -0.398, -0.448, -0.448, -0.090,
     + -0.078, -0.078, -0.040, -0.040,  0.025,  0.377,  0.377,
     +  0.152,  0.152, -1.316, -0.613, -0.496, -0.496, -0.628,
     + -0.783, -0.707, -0.707, -0.803, -0.596, -0.966, -0.937,
     + -0.680, -0.888, -0.869, -0.869, -0.974, -0.608, -0.675,
     + -0.655, -0.159, -0.396, -0.307, -0.307, -0.404, -0.218,
     + -0.218, -0.225,  0.432,  0.595,  0.595,  0.683,  0.857,
     +  0.400,  1.362,  1.362,  1.544,  1.843,  1.843,  1.548,
     +  1.662,  1.662,  1.467,  1.627,  1.419,  1.604,  1.325,
     +  1.381,  1.381,  0.762,  0.762,  1.372,  1.412,  1.633,
     +  1.472,  1.463,  1.527,  1.401,  1.710,  1.655,  1.655,
     +  1.004,  0.971,  1.245,  0.849,  0.847,  1.151,  0.876,
     +  1.047,  1.047,  1.209,  1.067,  1.067,  1.024,  1.024,
     +  0.973,  1.112,  1.164,  1.053,  1.053,  1.430,  1.430,
     +  1.329,  1.204,  1.480,  1.395,  1.543,  1.495,  1.260,
     +  1.260,  1.477,  1.477,  1.573,  1.519,  1.459,  1.308,
     +  1.308,  1.335,  0.964,  0.964,  1.084,  1.084,  1.192,
     +  0.545,  0.460,  0.462,  0.662,  0.662,  0.742,  0.647,
     +  0.647,  0.671,  0.751,  0.751,  0.588,  0.588,  0.421,
     +  0.421,  0.358,  0.358,  0.189,  0.189,  0.268,  0.268,
     +  0.183,  0.300,  0.300, -0.378, -0.462, -0.392, -0.392,
     + -0.395, -0.533, -0.693, -0.693, -0.676, -0.676, -0.364,
     + -0.351, -0.363, -0.363, -0.402, -0.382, -0.407, -0.407,
     + -0.292, -0.461, -0.461, -0.229, -0.229, -0.420, -0.474,
     + -0.474, -0.767, -0.501, -0.501, -0.869, -0.869, -0.504,
     + -0.708, -0.620, -0.620, -0.830, -0.651, -0.651, -0.856,
     + -0.624, -0.828, -0.828, -0.861, -0.861, -0.775, -0.716,
     + -0.716, -0.606, -0.798, -0.865, -0.867, -0.850, -0.850,
     + -0.893, -0.818, -0.846, -0.846,  0.213,  0.089,  0.227,
     +  0.085, -0.093, -0.180, -0.110, -0.573, -0.231, -0.005,
     + -0.431, -0.137, -0.137, -0.269, -0.269, -0.067, -0.067,
     + -0.070, -0.070,  0.028,  0.232, -0.257, -0.103, -0.020,
     +  0.026, -0.162, -0.010, -0.140, -0.009, -0.009, -0.096,
     + -0.107, -0.107, -0.191, -0.023, -0.305,  0.123,  0.009,
     + -0.168, -0.005,  0.170,  0.170, -0.078, -0.189, -0.189,
     + -0.034, -0.034,  0.129,  0.129, -0.034, -0.034, -0.205,
     + -0.263, -0.591, -0.591,  0.037, -0.071,  0.071,  0.071,
     + -0.175, -0.098, -0.139, -0.171, -0.352, -0.277, -0.277,
     + -0.417, -0.249, -0.285, -0.089,  0.415, -0.128, -0.128,
     + -0.216, -0.244, -0.349, -0.349, -0.289, -0.289, -0.299,
     + -0.340, -0.269, -0.293, -0.206, -0.213, -0.212, -0.030,
     + -0.124, -0.207, -0.173, -0.173, -0.043, -0.486, -0.486,
     + -0.075, -0.075, -0.053, -0.053, -0.303, -0.149,  0.000 /
      Data RunShf97_v1 /
     +  0.000, -0.799, -0.777, -0.777, -0.876, -0.797, -0.876,
     + -0.876, -1.303, -0.962, -0.962, -0.772, -0.576, -0.576,
     + -0.644, -0.644, -0.610, -0.622, -0.827, -0.718, -0.718,
     + -0.435,  0.004, -0.125, -0.137, -0.279, -0.341,  0.001,
     + -0.019, -0.019, -0.414, -0.414, -0.115, -0.300, -0.174,
     + -0.083, -0.083, -0.010, -0.150, -0.318, -0.433, -0.044,
     + -0.179, -0.390,  0.012, -0.201, -0.259, -0.214, -0.127,
     + -0.127, -0.226, -0.158, -0.321, -0.364, -0.356, -0.272,
     + -0.262, -0.262, -0.158, -0.158, -0.219, -0.219, -0.132,
     + -0.258, -0.241, -0.379, -0.092, -0.434, -0.434, -0.427,
     + -0.191, -0.087,  0.144, -0.332, -0.386, -0.186, -0.186,
     + -0.259, -0.322, -0.122, -0.229, -0.229, -0.306, -0.127,
     + -0.270, -0.254, -0.225, -0.221,  0.249, -0.169, -0.169,
     + -0.061,  0.106, -0.100,  0.181, -0.048,  0.034,  0.011,
     +  0.011,  0.248, -0.102, -0.102, -0.001,  0.031, -0.254,
     +  0.127,  0.125, -0.068,  0.008, -0.299, -0.224, -0.329,
     + -0.096, -0.030, -0.030, -0.206, -0.280, -0.124, -0.124,
     + -0.170, -0.118, -0.118, -0.166, -0.058, -0.058,  0.004,
     + -0.074,  0.052,  0.098,  0.098, -0.232,  0.088,  0.088,
     +  0.084, -0.024,  0.318, -0.110, -0.110, -0.075, -0.032,
     + -0.032, -0.325, -0.075,  0.066,  0.066, -0.282, -0.050,
     +  0.022, -0.185, -0.071, -0.071, -0.008,  0.014,  0.192,
     + -0.130, -0.130,  0.168, -0.244,  0.224,  0.090,  0.043,
     +  0.139, -0.075, -0.221, -0.221,  0.082,  0.143,  0.143,
     +  0.154,  0.154, -0.109,  0.046,  0.046,  0.024, -0.118,
     + -0.118,  0.064,  0.136,  0.062,  0.237,  0.116,  0.090,
     + -0.025,  0.301,  0.034,  0.034,  0.185, -0.058, -0.058,
     +  0.094,  0.071, -0.235, -0.084,  0.079, -0.145,  0.029,
     + -0.010, -0.104,  0.244,  0.244, -0.317, -0.317,  0.076,
     + -0.122, -0.004,  0.003, -0.031, -0.627, -0.019, -0.203,
     + -0.098,  0.055,  0.006, -0.015, -0.015, -0.026, -0.026,
     +  0.072,  0.099,  0.161,  0.161,  0.112, -0.091,  0.182,
     +  0.182, -0.014, -0.005,  0.088,  0.053, -0.141, -0.045,
     + -0.116,  0.492, -0.022, -0.015,  0.127,  0.118,  0.118,
     + -0.040, -0.040,  0.014,  0.243,  0.243,  0.178,  0.178,
     + -0.019,  0.024, -0.011,  0.035,  0.078,  0.014,  0.012,
     +  0.096,  0.077,  0.094,  0.109,  0.164,  0.255,  0.168,
     +  0.235,  0.530,  0.775,  0.746,  0.811,  0.337,  0.226,
     +  0.350,  0.252,  0.252,  0.348,  0.324,  1.315,  0.553,
     +  1.200,  1.210, -0.165,  0.509,  0.509,  0.380,  0.764,
     +  0.640,  0.528,  0.528,  0.630,  0.630,  0.788,  0.788,
     +  0.541,  0.662,  0.378,  0.378,  0.592,  0.592,  0.594,
     +  0.542,  0.542,  0.674,  0.409,  0.781,  0.745,  0.673,
     +  0.673,  0.769,  0.679,  0.751,  0.821,  0.911,  0.799,
     +  0.751,  0.689,  0.633,  0.738,  0.816,  0.912, -0.095,
     + -0.271, -0.135, -0.311, -0.195, -0.182, -0.145, -0.249,
     + -0.304, -0.304,  0.018,  0.184,  0.079,  0.058,  0.058,
     +  0.069,  0.069, -0.024, -0.090, -0.090,  0.168,  0.523,
     +  0.187,  0.145,  0.145,  0.159,  0.191,  0.318,  0.315,
     +  0.064, -0.016,  0.189,  0.244,  0.388,  0.391, -0.341 /
      Data RunShf97_v2 /
     +  0.752,  0.752,  0.742,  0.784,  0.640,  0.737,  0.439,
     +  0.535,  0.443,  0.993,  0.650,  0.274,  0.474,  0.474,
     +  0.251,  0.072, -0.478, -0.478, -0.362, -0.146, -0.146,
     + -0.527, -0.079,  0.135, -0.015, -0.242, -0.035, -0.035,
     +  0.106,  0.137,  0.137, -0.549, -0.585, -0.503,  0.061,
     +  0.232,  0.232, -0.047,  0.266,  0.266,  0.042, -0.015,
     + -0.001,  0.142,  0.247,  0.107,  0.030, -0.042, -0.239,
     + -0.239,  0.163,  0.484,  0.285,  0.010,  0.331,  0.331,
     + -0.055,  0.344,  0.285, -0.129,  0.209,  0.250,  0.250,
     +  0.192,  0.192,  0.088,  0.036,  0.185,  0.094, -0.006,
     +  0.002, -0.023, -0.023, -0.010, -0.137, -0.130, -0.200,
     + -0.259, -0.259, -0.109,  0.054,  0.054,  0.018, -0.149,
     + -0.264, -0.144, -0.010,  0.267,  0.156,  0.148,  0.231,
     +  0.211,  0.075,  0.657,  0.757,  0.160,  0.160,  1.063,
     +  0.774,  0.712,  0.876,  0.567,  0.556,  0.556,  0.420,
     +  0.578,  0.333,  0.243,  0.490,  0.480,  0.388,  0.357,
     +  0.258,  0.462,  0.462,  0.199,  0.199,  0.179,  0.334,
     +  0.175, -0.056,  0.129,  0.209,  0.209,  0.401,  0.401,
     +  0.444,  0.444,  0.571,  0.571,  0.295,  0.182,  0.416,
     +  0.206,  0.206,  0.362,  0.322,  0.145,  0.102, -0.019,
     +  0.260,  0.260,  0.228,  0.250,  0.012, -0.037, -0.222,
     +  0.200,  0.029, -0.040,  0.012,  0.166,  0.089,  0.089,
     + -0.011, -0.011, -0.175, -0.061, -0.160, -0.160,  0.063,
     +  0.018, -0.126,  0.250,  0.207,  0.077, -0.174, -0.101,
     +  0.273,  0.169,  0.230,  0.029,  0.169,  0.131,  0.131,
     + -0.345,  0.105,  0.105,  0.129,  0.110,  0.300,  0.054,
     +  0.054,  0.218, -0.072,  0.363,  0.147,  0.155,  0.398,
     +  0.295,  0.213,  0.117,  0.218,  0.224,  0.371,  0.363,
     +  0.384, -0.044,  0.326,  0.235, -0.087, -0.087,  0.209,
     +  0.026, -0.108, -0.229, -0.232, -0.232, -0.614, -0.614,
     + -0.547, -0.463, -0.535, -0.535, -0.496, -0.357, -0.357,
     + -0.773, -0.537, -0.537, -0.279, -0.323, -0.190, -0.578,
     + -0.210,  0.058,  0.058, -0.510, -0.493, -0.476, -0.292,
     + -0.086, -0.086, -0.006, -0.070, -0.113, -0.113, -0.156,
     + -0.156, -0.210, -0.176, -0.110, -0.110, -0.143, -0.143,
     + -0.354, -0.354, -0.320, -0.195, -0.239, -0.151, -0.301,
     + -0.477, -0.686, -0.648, -0.648, -0.244, -0.244,  0.855,
     + -0.261, -0.089, -0.106,  0.098,  0.098, -0.052, -0.052,
     + -0.099, -0.113, -1.575, -1.388, -1.426, -1.270, -1.294,
     + -1.433, -1.511, -1.385, -1.310, -1.068, -1.133, -1.133,
     + -1.013, -1.119, -1.171, -1.649, -1.155, -1.383, -1.458,
     + -0.306, -0.306, -0.268, -0.268, -0.232, -0.232, -1.306,
     + -1.306, -1.427, -1.232, -1.453, -1.313, -1.313, -1.108,
     + -1.189, -1.189, -1.201, -2.363, -2.363, -2.362, -2.406,
     + -2.642, -2.330, -2.330, -2.423, -2.423, -1.434, -0.573,
     + -0.573, -0.793, -0.929, -0.929, -0.977, -1.462, -1.408,
     + -1.416, -1.331, -1.428, -1.617, -1.742, -1.697, -1.746,
     + -1.673, -1.673, -1.729, -1.436, -1.598, -1.706, -1.706,
     + -1.872, -1.872, -1.813, -1.765, -1.632, -1.650, -1.887,
     + -1.801, -1.801, -1.518, -1.785, -1.933, -1.932, -1.632 /
      Data RunShf97_v3 /
     + -1.808, -1.638, -1.505, -1.582, -1.473, -1.299, -1.299,
     + -1.079, -1.406, -1.406, -1.073, -1.073, -0.994, -0.994,
     + -0.949, -1.001, -0.984, -1.009, -0.839, -0.874, -0.850,
     + -0.850, -1.235, -1.235, -0.973, -0.979, -1.216, -1.216,
     + -0.953, -0.984, -1.192, -0.697, -0.697, -0.738, -0.883,
     + -0.883, -0.936, -0.655, -0.683, -0.846, -0.846, -0.699,
     + -0.548, -0.731, -0.597, -0.613, -0.922, -0.922, -0.678,
     + -0.904, -0.904, -0.877, -0.806, -0.865, -0.895, -0.943,
     + -0.975, -0.766, -0.766, -0.825, -0.900, -1.456, -1.323,
     + -1.780, -1.780, -1.710, -1.710, -2.025, -2.025, -1.748,
     + -1.748, -1.066, -1.573, -1.302, -1.690, -1.384, -1.547,
     + -1.692, -1.508, -1.508, -1.646, -1.646, -1.761, -1.534,
     + -1.534, -1.615, -1.593, -1.702, -1.731, -1.646, -1.317,
     + -1.568, -1.568, -1.405, -1.606, -1.606, -1.693, -1.797,
     + -1.885, -1.613, -1.708, -1.618, -1.845, -1.650, -1.650,
     + -1.589, -1.589, -1.998, -1.897, -1.803, -1.603, -1.652,
     + -1.652, -1.610, -1.585, -1.448, -1.790, -1.822, -1.720,
     + -1.685, -1.671, -1.728, -1.719, -1.741, -1.826, -1.826,
     + -1.576, -1.576, -1.557, -1.557, -1.008, -0.958, -0.623,
     + -0.821, -0.811, -0.811, -0.813, -0.644, -0.732, -0.979,
     + -0.912, -0.857, -1.110, -1.110, -0.692, -0.600, -0.533,
     + -0.534, -0.402, -0.608, -0.581, -0.581, -0.752, -0.752,
     + -0.820, -0.820, -0.523, -0.523, -0.725, -0.782, -0.620,
     + -0.397, -0.554, -0.537, -0.585, -0.841, -0.749, -0.749,
     + -1.905, -2.075, -2.037, -2.004, -2.004, -1.920, -1.920,
     + -0.532, -0.617, -0.522, -0.597, -0.556, -0.659, -0.460,
     + -0.460, -0.387, -0.702, -0.702, -0.544, -0.506, -0.386,
     + -0.792, -0.724, -0.557, -0.557, -0.354, -0.585, -0.596,
     + -0.498, -0.622, -0.809, -0.553, -0.666, -0.571, -0.656,
     + -0.648, -0.648, -0.622, -1.012, -1.012, -0.686, -0.686,
     + -0.786, -0.753, -0.725, -0.667, -0.713, -0.512, -0.657,
     + -0.657,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000 /

C
      If(AccRun.eq.run) Then
         TfRunCorr = AccShf
         GoTo 222
      EndIF
      AccRun = run
C 97
      If(run.gt.25000) Then
          Do i=1,MaxRCorr97
             If(run.le.RunNum97(i)) Then
               inum=i
               GoTo 101
             EndIf
          EndDo
 101      Continue
          val1=RunShf97(inum-1)
          run1=float(RunNum97(inum-1))
          val2=RunShf97(inum)
          run2=float(RunNum97(inum))
C 96
      ElseIf(run.gt.20000) Then
          Do i=1,MaxRCorr96
             If(run.le.RunNum96(i)) Then
               inum=i
               GoTo 102
             EndIf
          EndDo
 102      Continue
          val1=RunShf96(inum-1)
          run1=float(RunNum96(inum-1))
          val2=RunShf96(inum)
          run2=float(RunNum96(inum))
C 95
      ElseIf(run.gt.10330) Then
          Do i=1,MaxRCorr95
             If(run.le.RunNum95(i)) Then
               inum=i
               GoTo 103
             EndIf
          EndDo
 103      Continue
          val1=RunShf95(inum-1)
          run1=float(RunNum95(inum-1))
          val2=RunShf95(inum)
          run2=float(RunNum95(inum))
C 94
      Else
          Do i=1,MaxRCorr94
             If(run.le.RunNum94(i)) Then
               inum=i
               GoTo 104
             EndIf
          EndDo
 104      Continue
          val1=RunShf94(inum-1)
          run1=float(RunNum94(inum-1))
          val2=RunShf94(inum)
          run2=float(RunNum94(inum))
      EndIf
C
      rrun=float(run)
      TfRunCorr = ((rrun-run1)/(run2-run1))*(val2-val1)
     +           + val1
      AccShf = TfRunCorr
C
C      Write(*,*) ' run - to - run correction = ',AccShf
C
 222  Continue
C
      Return
      End
