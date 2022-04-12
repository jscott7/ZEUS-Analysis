C     =============================================
      Subroutine ZufosNT( Vtx,  NrEcells, CellList,
     +                    ElPos,ElProb,   ElEnergy)
C     =============================================
c     INPUT: Vtx(1:3)  :
c            NrEcells  : nr of cells from electron
c            Celllist  : Array with cellnumbers
c            Elpos(1:3): x,y,z position of electron 
c            Elprob    : electron probability
c            ElEnergy  : electro energy
c
c     This routine is an overhead routine, which initialises
c     the zufos properly, gets the timing information for calorimeter
c     islands, finds backsplash, and in the end loops over all zufos
c     to add Px,Py,Pz,E,...
c     It fills the Common block zRecNT01;
c
c     Px,   Py,   Pz,   E:     sum of all zufos
c     PxCal,PyCal,PzCal,ECal:  sum of the calorimeter part of the zufos
c                              - 1st index: F/R/BCAL (1/2/3)   
c                              - 2nd index: EMC/HAC  (1/2)
c     PxTrk,PyTrk,PzTrk,ETrk:  sum of tracking part of the zufos
c     PxBSp,PyBSp,PzBSp,EBSp:  sum of all BackSplash zufos;
c                              NOT YET SUBTRACTED FROM THE TOTAL SUM
c                              - index: 3 different criteria for
c                                gammamax. 
c                                2 SHOULD BE USED,
c                                1&3 ARE MENT FOR SYSTEMATIC CHECKS 
c     EminPzNT,PtNT,GammaNT: Ready-to-use-Take-it-or-leave-it-quantities
c
c     v8.0: Bug Fix in Calorimeter Islanding!
c     ------------------------------------------------------------------
c
c     ----------------------------
c     Author: N.Tuning 10/2/99
c     Zufo Version: 4.6  (ZufosNT)
c     17/06/99  Change wrt version 4: 
c               4.6 --> 6.2 Corrected : BUG
c     1) The timing of the Calorimeter Islands was calculated
c        with one PMT side only.
c     2) Arrays were not checked for their length. This causes
c        problems only if user calls ZufosNT with very large Zvtx.
c     3) Included check on Zvtx 
c
c     --------------------------------
c     23/07/00  Change wrt version 6: 
c               6.2 --> 7.0
c     1) Include ET variables
c     2) Give maximum size of array's in AvgTimClu (not via argument!)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit None
c
#include "partap.inc"
#include "caltru.inc"
#include "zdskey.inc"
c
#include "zisles.inc"
#include "zdrecgb.inc"
#include "ZufosNT.inc"
c
      Integer NrEcells, CellList(*)
      Integer Ierr, Icell, I, J
      Integer MoreTr, SecTr, EneCor, BSpCor
      Real    ElPos(3), ElProb,   ElEnergy
      Real    Vtx(3), xv(3)
      Real    Pt2
      Integer Iwarning
      Logical First
      Data    First/.TRUE./

      If (First) Then
         Iwarning = 0.
         First    = .FALSE.

         CALL Z_INI_ISLES

         Write(*,*) ' -------------------------------------------------'
         Write(*,*) ' This is the improved version of Genadys zufos    '
         Write(*,*) '                                                  '
         Write(*,*) '          ZufosNT  Version v8.0  22/05/00         '
         Write(*,*) '          ===============================         '
         Write(*,*) '                                                  '
         Write(*,*) ' Documentation:                                   '
         Write(*,*) ' www-zeus.desy.de/~ntuning/ZEUS_ONLY/zufo_new.html'
         Write(*,*) ' -------------------------------------------------'

      EndIf

      xv(1) = Vtx(1)
      xv(2) = Vtx(2)
      xv(3) = Vtx(3)

      If ( Abs( Vtx(1)).GT.10. .OR. Abs( Vtx(2)).GT.10.) Then
         Write(*,*) '*** WARNING: x or y vertex bigger then 10 cm.'
         Write(*,*) '             Check x,y and z vertex in ZufosNT!'
         Write(*,*) '             --> Put vtx(3) at (0,0,0)'
         xv(1) = 0.
         xv(2) = 0.
         xv(3) = 0.
      EndIf

      If ( Vtx(3).LT.-148. .OR.  Vtx(3).GT.235. ) Then
         Write(*,*) '*** WARNING: z vertex out of range (-148.:235.)'
         Write(*,*) '   (ZufosNT) --> Put vtx(3) at (0,0,0)'
         xv(1) = 0.
         xv(2) = 0.
         xv(3) = 0.
      EndIf



C ----------------------
C --  Initialise output:
C ----------------------
      Px = 0.
      Py = 0.
      Pz = 0.
      E  = 0.
      Et = 0.
      PxTrk = 0.
      PyTrk = 0.
      PzTrk = 0. 
      ETrk  = 0.
      EtTrk = 0.
      Call Vzero( PxCal, 6)
      Call Vzero( PyCal, 6)
      Call Vzero( PzCal, 6)
      Call Vzero(  ECal, 6)
      Call Vzero( EtCal, 6)
      Call Vzero( PxBSp, 3)
      Call Vzero( PyBSp, 3)
      Call Vzero( PzBSp, 3)
      Call Vzero(  EBSp, 3)
      Call Vzero( EtBSp, 3)

      EminPzNT = 0.
      Pt2      = 0.
      PtNT     = 0.
      GammaNT  = 0.
      EtNT     = 0.

C --------------------------------
C --- Fill zRec02 COMMON block ---
C --------------------------------
      eNdis  = NrEcells
      ePrdis = ElProb
      eEdis  = ElEnergy
      eXdis  = ElPos(1)
      eYdis  = ElPos(2)
      eZdis  = ElPos(3)
 
C --------------------------------
C --- Fill zIsl01 COMMON block ---
C --------------------------------
      zCells = CouTab(CalTru)
      If (zCells.GT.Nr_cMax) Then
         CALL FETTAB(ZDSKEY,ID,1)
         Write(6,*)'*==> ',ZDSKEY_NR1, ZDSKEY_NR2,
     &        'Caltru has more than',Nr_cMax,
     &        ' cells:',zCells,' cell in this event'
         RETURN
      EndIf

      Icell = 0
      Do I=1,zCells
         Call FetTab(Caltru,ID,I)
         Do J=1,MIN(eNdis,50)
            If (Caltru_Cellnr.EQ.CellList(J)) GOTO 110
         EndDo
         Icell        = Icell + 1
         zPnrl(Icell) = Caltru_CellNr
         zEl(  Icell) = Caltru_e
         zImbl(Icell) = Caltru_imbal
         zID(  Icell) = Caltru_ID
         zT1(  Icell) = Caltru_t(1)
         zT2(  Icell) = Caltru_t(2)
 110     Continue
      EndDo
      zCells = Icell
 

C -----------------
C --  Initialise:
C -----------------
      MoreTr = 1 
      SecTr  = 1 
      EneCor = 0 
      BSpCor = 1
      Call z_Ini_RecGB( MoreTr, SecTr, EneCor, BSpCor)
 
C --------------------------
C --  Do the ZUFO thing:
C --------------------------
      Call z_RecGB( 1,1,xv,Ierr)

C -------------------------------------------------------------
C --  Get the timing information for the calorimeter islands:
C --------------------------------------------------------------
      Call Isles_Time

C ----------------------------------------------------------
C --  Flag the zufos which are considered to be BackSplash:
C ----------------------------------------------------------
      Call BcorZufo

C --------------------------------------------------------
C --  Loop over all zufos and make sums of Px,Py,Pz,E :
C --------------------------------------------------------
      Call ZufoLoop

      EminPzNT = E - Pz - (EBSp(2)-PzBSp(2))
      Pt2      = (Px-PxBSp(2))**2 + (Py-PyBSp(2))**2
      PtNT     = Sqrt( Pt2)
      EtNT     = Et - EtBSp(2)

      IF (EminPzNT.GT.0. .AND. Pt2.GT.0.) THEN
         GammaNT  = ACOS( (Pt2 - EminPzNT**2)/(Pt2 + EminPzNT**2))

      ELSEIF (Iwarning.LE.10) THEN
         GammaNT  = 0.
         Write(*,*) '*** WARNING (ZufosNT): E-Pz(had), Pt(had):',
     +        EminPzNT, PtNT
         Iwarning = Iwarning + 1
         If (Iwarning.EQ.10) Write(*,*) 
     +        'Warning (ZufosNT) occured 10 times, no print anymore...'

      ELSE
         GammaNT  = 0.

      ENDIF

      END

C     -----------------------------------------------------------------
C     -----------------------------------------------------------------

C     ========================
      Subroutine Isles_Time( )
C     ========================
c     Calls:   AvgTimClu
c     Purpose: Get timing information for islands in "zisles.inc".
c     Author:  Niels Tuning 11/98
c              
c     N.T. Modified 11/2/99: added some boundary checks.
c     N.T. Modified 12/06/99: corrected bug! Tcell2 was not filled!!
c
      Implicit None

#include "zisles.inc"
#include "isles_time.inc"

      Integer NrCellIsl(nI_Max)
      Integer I
      Integer Icell, kIsl


      Do kIsl=1,nI_Max
         tIsl(  kIsl)    = 0.
         teIsl( kIsl)    = 0.
         NrCellIsl(kIsl) = 0
      Enddo
     
      IF (nIsl .GT. nI_Max) THEN
         WRITE(*,*) 
     +        'ERROR Isles_Time: Nr of islands too big (', nI_Max,') :',
     +        nIsl
         RETURN
      ENDIF
      IF (zCells .GT. Nr_cMax) THEN
         WRITE(*,*) 
     +        'ERROR Isles_Time: Nr of cells too big (',Nr_cMax,')  :',
     +        zCells
         RETURN
      ENDIF

      Do Icell=1,zCells
         kIsl = I_caltru(Icell)
         IF (kIsl .GT. 0) THEN
            NrCellIsl(                 kIsl) = NrCellIsl(kIsl) + 1
            TCellIsl1( NrCellIsl(kIsl),kIsl) = zT1(  Icell) 
            TCellIsl2( NrCellIsl(kIsl),kIsl) = zT2(  Icell) 
         ELSE
            WRITE(*,*) 
     +           'ERROR Isles_Time: Caltru ID of cell ',Icell,
     +           '  is zero ? :  ',I_caltru(Icell),'   return...'
            RETURN
         ENDIF   
      EndDo
      
      Do kIsl=1,MAX( nIsl, nI_Max)
         NrCells = NrCellIsl( kIsl)

         Do I=1,MAX( NrCellIsl( kIsl),Nr_cIslMax)
            TCell1( I) = TCellIsl1( I,kIsl)
            TCell2( I) = TCellIsl2( I,kIsl)
            ECell(  I) = ECellIsl(  I,kIsl)
            ImbCell(I) = ImbCellIsl(I,kIsl)
         EndDo
         Call AvgTimClu(NrCells,TCell1,TCell2,ECell,ImbCell,Time,Timerr)
         tIsl( kIsl) = Time
         teIsl(kIsl) = Timerr
      EndDo
      End


C     ==============================================================
      Subroutine AvgTimClu(NrCells,TCell1,TCell2,ECell,ImbCell,
     +                     Time,Timerr)
C     ==============================================================
c     Purpose: Get the average timing for NrCells. 
c     Author:  Niels Tuning 11/98
c
C     v7.0 23/03/00 N.Tuning Change declaration of array.     
C
      Implicit None
      
      Integer   Nr_cMax
      Parameter( Nr_cMax=1000)
c
      Integer   I
      Integer   NrCells
      Real      TCell1( Nr_cMax), TCell2( Nr_cMax)
      Real      ECell(  Nr_cMax)
      Real      ImbCell(Nr_cMax)
      Real      Time_tmp, Time, Timerr
      Real      Epmt1, Epmt2, twgt1, twgt2, twtot

      Time_tmp = 0.
      twtot    = 0.

      Do I=1,NrCells
         twgt1 = 0.
         twgt2 = 0.
         Epmt1 = 0.
         Epmt2 = 0.

         If (ImbCell(I).eq.0.) Then
            Epmt1 = 0.5*ECell(I)
            If (Epmt1.gt.0.1) twgt1 = 1./(0.4 + 1.4/Epmt1**0.65)**2
            twtot = twtot + twgt1
            Time_tmp = Time_tmp + TCell1(I)*twgt1
         Else
            Epmt1 = (ECell(I) - ImbCell(I))/2
            Epmt2 = (ECell(I) + ImbCell(I))/2
            If (Epmt1.gt.0.1) twgt1 = 1./(0.4 + 1.4/Epmt1**0.65)**2
            If (Epmt2.gt.0.1) twgt2 = 1./(0.4 + 1.4/Epmt2**0.65)**2
            twtot = twtot + twgt1 + twgt2 
            Time_tmp = Time_tmp + TCell1(I)*twgt1
            Time_tmp = Time_tmp + TCell2(I)*twgt2
         EndIf

      EndDo

      If (twtot.gt.0.) Then
         Time   = Time_tmp/twtot
         Timerr = 1./sqrt( twtot)
      Else
         Time   = 0.
         Timerr = 0.
      EndIf

      End

C     -----------------------------------------------------------------
C     -----------------------------------------------------------------

C     =====================
      Subroutine Bcorzufo()
C     =====================
c     This routine loops over all zufos and flags the zufos which
c     are considered to be backsplash.
c
c     o Only calorimeter islands without tracks pointing to it
c       (zufo typ 31) are considered.
c     o Below 3 GeV
c     o Islands in the FCAL are never considered backsplash.
c     o Islands with very good timing (E>1.5 GeV , abs(time)<5 ns)
c       are not considered backsplash.
c     o Islands behind  gammamax=gamma+50 deg are considered as 
c       backsplash.
c
c     The procedure is iterative. After each iteration, gamma
c     changes and thus the criteria to cut backsplash.
c     THIS STOPS WHEN GAMMA CHANGES LESS THEN 1%,
c     but never more then 3 iterations.      
c
c     The biggest effect of all the criteria above, is the 
c     functional form of gammamax. This is the reason that this 
c     is varied in order to do systematic studies.
c
c     Author: N.Tuning 22/1/99
c     Zufo Version: 6.2  (ZufosNT) 
c     Note that the timing is now correct!!
c     -------------------------------------------------------
      Implicit None
c
#include "zdrecgb.inc"
#include "zisles.inc"
#include "isles_time.inc"
#include "ZufosNT.inc"
c
      Integer I, ZufoTyp, TrackId, Niter, VaryCrit
      Real    CalEner, CalTime, CalTimerr, Theta, CalZpos
      Real    Px_tmp, Py_tmp, Pz_tmp, E_tmp, Pt2, Epz
      Real    gamma(4), gamma0, gammac, gammamax
      Real       Emax
      Parameter (Emax = 3.)

      Logical Take
      Integer CountErr
      Logical First
      Data    First/.TRUE./

      If (First) Then
         First = .FALSE.
         CountErr = 0
      EndIf

      Do i=1,zufos_MaxNT
         BSpFlag(1,I) = 0
         BSpFlag(2,I) = 0
         BSpFlag(3,I) = 0
      EndDo

C     ============================
      IF (.NOT. BcorNT) GOTO 2000
C     ============================

      Px_tmp        = 0.
      Py_tmp        = 0.
      Pz_tmp        = 0.
      E_tmp         = 0.
      Do I=1,Nzufos
         Px_tmp = Px_tmp + zufo(1,I)
         Py_tmp = Py_tmp + zufo(2,I)
         Pz_tmp = Pz_tmp + zufo(3,I)
         E_tmp  = E_tmp  + zufo(4,I)
      ENDDO
      Pt2    = Px_tmp**2 + Py_tmp**2
      Epz    = E_tmp - Pz_tmp

      If (Pt2+Epz**2 .GT. 0.) Then
         gamma0 = ACOS( (Pt2-Epz**2)/(Pt2+Epz**2) )
      Else
         CountErr = CountErr + 1
c         If (CountErr.le.10) 
c     +        Write(*,*) '*** Bcorzufo; tering; no hadr. activity',
c     +                    CountErr, Pt2+Epz**2
         RETURN
      EndIf

C     -----------------------------------------------------------------
      
      Do VaryCrit = 1,3

         Niter = 0
 1000    Continue
         Niter = Niter + 1
         
            
         If (Niter.eq.1) Then
            gammac=gamma0
         ElseIf (Niter.gt.1) Then
            gammac=gamma(Niter-1)
         EndIf

         If (VaryCrit.eq.1) gammamax = gammac + 0.873 - 0.2
         If (VaryCrit.eq.2) gammamax = gammac + 0.873 
         If (VaryCrit.eq.3) gammamax = gammac + 0.873 + 0.2

         Px_tmp = 0.
         Py_tmp = 0.
         Pz_tmp = 0.
         E_tmp  = 0.
         Do I=1,Nzufos
            Take         = .TRUE.
            ZufoTyp      = tufo(1,I)
            Theta        = atan2(sqrt(zufo(1,I)**2 + zufo(2,I)**2),
     +                                zufo(3,I)                    )
            If (tufo(2,I).gt.0) Then
               TrackId   = vcthID(tufo(2,I))
            Else
               TrackId   = 0
            EndIf
            If (tufo(3,I).gt.0) Then
               CalTime   = tIsl( tufo(3,I))
               CalTimerr = teIsl(tufo(3,I))
               CalZpos   = zIsl( tufo(3,I))
            Else
               CalTime   = -999.
               CalTimerr = -999.
               CalZpos   = -999.
            EndIf
            
C--   ZUFO from Tracking?
            If ( ZufoTyp .LT. 30 .AND.
     +           ZufoTyp .GT. 31      ) Then
C              zufo from secundary vtx track?
               If ( TrackId .LT. 0) Then
                  Take = .FALSE.
               EndIf
               
C--   ZUFO from Cal with some track match?
            ElseIf (ZufoTyp .EQ. 30) Then
               
C--   ZUFO from Cal without track match?
            ElseIf (ZufoTyp .EQ. 31) Then
               CalEner = zufo(4,I)
               If (CalEner.LT.Emax .AND. Theta.GT.gammamax) Then
                  Take = .FALSE.
               EndIf
C              save FCAL islands: 
               If (CalZpos.GT.220.)                Take = .TRUE.
C              save islands with >1.5 GeV and good timing:
               If ( CalEner   .GT.   1.5   .AND. 
     +              abs(CalTime).LT. 5.E-9 .AND. 
     +              CalTimerr .LT.   2.    .AND. 
     +              CalTimerr .GT.   0.          ) Take = .TRUE.
            EndIf
            
            If (Take) Then
               Px_tmp = Px_tmp + zufo(1,I)
               Py_tmp = Py_tmp + zufo(2,I)
               Pz_tmp = Pz_tmp + zufo(3,I)
               E_tmp  = E_tmp  + zufo(4,I)
            ElseIf (.NOT.Take .AND. BcorNT) Then
               BSpFlag( VaryCrit, I) = 1
            EndIf

         ENDDO         
C--      end loop over all zufos.
         
         Pt2          = Px_tmp**2 + Py_tmp**2
         Epz          = E_tmp - Pz_tmp
         gamma(Niter) = ACOS( (Pt2-Epz**2)/(Pt2+Epz**2) )
         
         If (Niter .eq. 1) Then
            IF (  ABS(gamma0-gamma(1))/gamma(1).GE.0.01 ) THEN
               GOTO 1000
            ENDIF
         ElseIf (Niter .gt. 1) Then
            IF (  ABS(gamma(Niter-1)-gamma(Niter))/gamma(Niter).GE.0.01 
     +           .AND.    Niter.LE.3                             ) THEN
               GOTO 1000
            ENDIF
         EndIf


      EndDo
C---  end loop over different gammamax criteria.

 2000 CONTINUE
      END

C     -----------------------------------------------------------------
C     -----------------------------------------------------------------

C     =======================================================
      Subroutine Ecorzufo(MC,ZufoTyp,CalEner,EMCEner, Factor)
C     =======================================================
c    Input: MC      logical, true if mc, false if data
c           ZufoTyp integer, tufo(1,nzufos)
c           CalEner real,    calorimeter energy in island 
c           EMCEner real,    energy in EMC section
c
c    Output:Factor real, correction factor for zufo(i,Nzufos)
c
c     Called by z_RecGB
c
c     Cal and emc energy are obtained from Common/zIsl02/ 
c     in file zisles.inc.
c     The IslandId for a given zufo is given by the variable
c     tufo(3,Nzufos)
c     E.g. EMCEner = emcEIsl( tufo(3,Nzufos))
c
c     Separate corrections are applied on EMC only, HAC only,
c     or mixed islands. 
c     No correction is applied on islands below 1.5 GeV
c     No corection is applied on the tracking; only     
c     typ>29 zufo's are corrected.    
c
c     Author: N.Tuning 23/12/98
c
c     23/03/00 REMINDER: 
c         Default: THIS ENERGY CORRECTION IS NOT APPLIED
c         Advice:  DO NOT USE IT ... 
c     -------------------------------------------------------
      Implicit None
c
      Logical MC
      Integer ZufoTyp
      Real    CalEner, EMCEner, HACEner, Factor
      Integer EmcHac
      Real    mipcor
c
      COMMON/NIELSPAR/FITPAR(4)
      REAL FITPAR
      REAL CORFUN
c
      REAL FITPARDA(4,3)
      REAL FITPARMC(4,3)
c
      DATA FITPARMC /   0.0000,  -1.0522,  -0.2969,   0.1977,
     &                  0.4885,  -3.5162,  -0.3148,   0.2038,
     &                  0.8775,   0.2884,   3.0650,   0.1000/
      DATA FITPARDA /   0.0758,  -1.1981,  -0.2645,   0.2109,
     &                  2.9402,  -2.6168,  -0.2760,   0.2020,
     &                  0.9558,   0.0987,   2.8537,   0.1000/
C     -------------------------------------------------------
      Factor  = 1.
      HACEner = CalEner - EMCEner

      If (ZufoTyp.LT.30 ) GOTO 999
      If (CalEner.GT.1.5) GOTO 999

      IF (    HACEner.lt.0.05 .AND. EMCEner.ge.0.05) THEN
         EmcHac=1
      ELSEIF (HACEner.ge.0.05 .AND. EMCEner.lt.0.05) THEN
         EmcHac=2
      ELSEIF (HACEner.ge.0.05 .AND. EMCEner.ge.0.05) THEN
         EmcHac=3
      ELSEIF (HACEner.lt.0.05 .AND. EMCEner.lt.0.05) THEN
         GOTO 999
      ENDIF

      If (EmcHac.eq.1 .and. ZufoTyp.eq.30) Then
c there are tracks, so not a photon! Apply mixed correction!
         EmcHac=3
      EndIf

      IF (MC) THEN
         FITPAR(1)=fitparmc(1,EmcHac)
         FITPAR(2)=fitparmc(2,EmcHac)
         FITPAR(3)=fitparmc(3,EmcHac)
         FITPAR(4)=fitparmc(4,EmcHac)
      ELSE
         FITPAR(1)=fitparda(1,EmcHac)
         FITPAR(2)=fitparda(2,EmcHac)
         FITPAR(3)=fitparda(3,EmcHac)
         FITPAR(4)=fitparda(4,EmcHac)
      ENDIF
C
C
      IF (EmcHac.eq.1) THEN       ! EMC
         Factor = CORFUN(log10(CalEner))
         mipcor = 0.
         If ( log10(CalEner) .lt. 0.) Then
            mipcor = 0.1+0.1*log10(CalEner)
         EndIf
         Factor = ((CalEner-mipcor)/CalEner)*Factor
      ELSEIF (EmcHac.eq.2) THEN   ! HAC
c        HAC energy only?? EMC below threshold -> add some...
         Factor = CORFUN(log10(CalEner))

c changed 22/1/99 NT.; the following gives a big overcorrection... v3.0
c         Factor = Factor*(CalEner+0.15)/CalEner 

      ELSEIF (EmcHac.eq.3) Then ! MIX
         Factor = FITPAR(1)
     +        +FITPAR(2)*log10(CalEner)
     +        +FITPAR(3)*(log10(CalEner))**2
      ENDIF   

 999  CONTINUE
      END
C    ====================
      REAL FUNCTION CORFUN(X)
C    ====================
      IMPLICIT NONE
      COMMON/NIELSPAR/FITPAR(4)
      REAL FITPAR
      REAL X, FIX
      IF (X.GE.FITPAR(3)) THEN
         CORFUN = 1.+FITPAR(4)*EXP(-1.*FITPAR(1)*X)
      ELSEIF (X.LT.FITPAR(3)) THEN
         FIX   = 1. + FITPAR(4)*EXP(-1.*FITPAR(1)*FITPAR(3))
     +              - FITPAR(2)*FITPAR(3)
         CORFUN = FIX+FITPAR(2)*X
      ENDIF
      END

C     -----------------------------------------------------------------
C     -----------------------------------------------------------------

C     ===================
      Subroutine ZufoLoop
C     ===================
c
c     Author: N.Tuning 11/1/99
c     -------------------------------------------------------
      Implicit None
c
#include "zdrecgb.inc"
#include "zisles.inc"
#include "ZufosNT.inc"
#include "partap.inc"
#include "fmckin.inc"
c
      Integer I, IslId, IFRB, VaryCrit
      Real    EmcFrac, HacFrac
c
c --- Energy Correction
      Logical MC
      Integer ZufoTyp
      Real    CalEner, EmcEner, Factor
c
      INTEGER   Ierr, vector(5)
      CHARACTER KIND*5

      DO I = 1, Nzufos

         Factor = 1.
         If ( EcorNT .AND. tufo(3,I).gt.0) Then
            MC      = (coutab(FMCKIN).gt.0)
            ZufoTyp = tufo(1,I)
            CalEner = eIsl(    tufo(3,I))
            EMCEner = emcEIsl( tufo(3,I))
            Call Ecorzufo(MC, ZufoTyp, CalEner, EMCEner, Factor)
         EndIf

C-- Add total 4-vector of event:
         Px = Px + Factor*zufo(1,I)
         Py = Py + Factor*zufo(2,I)
         Pz = Pz + Factor*zufo(3,I)
         E  = E  + Factor*zufo(4,I)
         Et = Et + SQRT( (Factor*zufo(1,I) )**2 + 
     +                   (Factor*zufo(2,I) )**2  )

C-- Add 4-vector in the calorimeter sections:
         If ( tufo(1,I).GE.30 ) Then
            IslId = tufo(3,I)
            If (     zI1(IslId) .GE. 220.) Then
               IFRB = 1
            ElseIf ( zI1(IslId) .LE. -140.) Then
               IFRB = 2
            Else
               IFRB = 3
            EndIf
            EmcFrac = Factor*(             emcEIsl( IslId))/eIsl( IslId)
            HacFrac = Factor*(eIsl( IslId)-emcEIsl( IslId))/eIsl( IslId)
            
            If (     tufo(1,I).EQ.30 .OR. tufo(1,I).EQ.31) Then
               PxCal(IFRB,1) = PxCal(IFRB,1) + EmcFrac*zufo(1,I)
               PyCal(IFRB,1) = PyCal(IFRB,1) + EmcFrac*zufo(2,I)
               PzCal(IFRB,1) = PzCal(IFRB,1) + EmcFrac*zufo(3,I)
               ECal( IFRB,1) = ECal( IFRB,1) + EmcFrac*zufo(4,I)
               EtCal(IFRB,1) = EtCal(IFRB,1) + 
     +              SQRT( (EmcFrac*zufo(1,I))**2+(EmcFrac*zufo(2,I))**2)
               PxCal(IFRB,2) = PxCal(IFRB,2) + HacFrac*zufo(1,I)
               PyCal(IFRB,2) = PyCal(IFRB,2) + HacFrac*zufo(2,I)
               PzCal(IFRB,2) = PzCal(IFRB,2) + HacFrac*zufo(3,I)
               ECal( IFRB,2) = ECal( IFRB,2) + HacFrac*zufo(4,I)
               EtCal(IFRB,2) = EtCal(IFRB,2) + 
     +              SQRT( (HacFrac*zufo(1,I))**2+(HacFrac*zufo(2,I))**2)
            ElseIf ( tufo(1,I).EQ.37 .OR. tufo(1,I).EQ.41) Then
               ECal( IFRB,1) = ECal( IFRB,1) + EmcFrac*zufo(4,I)
               ECal( IFRB,2) = ECal( IFRB,2) + HacFrac*zufo(4,I)
               PxTrk         = PxTrk + zufo(1,I)
               PyTrk         = PyTrk + zufo(2,I)
               PzTrk         = PzTrk + zufo(3,I)
               EtTrk         = EtTrk + SQRT( zufo(1,I)**2+zufo(2,I)**2)
            EndIf
         Else
            PxTrk         = PxTrk + zufo(1,I)
            PyTrk         = PyTrk + zufo(2,I)
            PzTrk         = PzTrk + zufo(3,I)
            ETrk          = ETrk  + zufo(4,I)
            EtTrk         = EtTrk + SQRT( zufo(1,I)**2+zufo(2,I)**2)
         EndIf

C-- Add 4-vector of the backsplash islands:
         Do VaryCrit = 1, 3
            If ( BSpFlag(VaryCrit,I).eq.1) Then
               PxBSp(VaryCrit) = PxBSp(VaryCrit) + Factor*zufo(1,I)
               PyBSp(VaryCrit) = PyBSp(VaryCrit) + Factor*zufo(2,I)
               PzBSp(VaryCrit) = PzBSp(VaryCrit) + Factor*zufo(3,I)
               EBSp( VaryCrit) = EBSp( VaryCrit) + Factor*zufo(4,I)
               EtBSp(VaryCrit) = EtBSp(VaryCrit) + 
     +              SQRT( (Factor*zufo(1,I)**2) + (Factor*zufo(2,I)**2))
            EndIf
         EndDo

      ENDDO
C     end loop over all zufos

      END
C     -----------------------------------------------------------------
