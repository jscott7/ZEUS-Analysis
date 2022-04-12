      Subroutine z_CAL_and_CTD(Zvtx,Ierr)
      Implicit NONE
*-------------------------------------------------------------------------------
*
*   Called by z_RecGB routine to build CAL_and_CTD ZUFOs. 
*
*     INPUT:  
*
*         Zvtx(3)  <--> vertex(1=x,2=y,3=z)
*
*     OUTPUT: 
*       
*       Ierr:   0  <--> No problemos
*
*   Author:   Gennady Briskin, Tel Aviv University
*   Date:     20-Dec-1996
*
*_______________________________________________________________________________
*
#include "zisles.inc"
#include "zdrecgb.inc"

*==>  Input
      Real       Zvtx(3)

*==>  OutPut
      Integer    Ierr

*==>  Useful Constants
      REAL       Pi, TwoPi
      PARAMETER (Pi = 3.14159265358979323846)
      PARAMETER (TwoPi = 2.*Pi)

      Real       piMass
      Parameter (piMass = 0.13956995)
      Real       xpi2
      Parameter (xpi2   = piMass**2)

*==>  Track-Island matching(VCTDCA)
      Real    DistCA,DCATerr,DCACerr
      Real    Coord(3),DCoord(3)
      Real    thIsl,hEIsl,epsIsl
      Real    matPr

      Real    Tr_E,E_Isl
      Real    DCA_hCut

*==>  Island Info
      Integer  hnI         , hcI(   nI_Max)
      Real     heI( nI_Max), emcI(  nI_Max)
      Real     hxI( nI_Max), hyI(   nI_Max)
      Real     hzI( nI_Max), hrI(   nI_Max)
      Real     hthI(nI_Max), hphI(  nI_Max)
      Integer  ITRm(nI_Max),ITrw(15,nI_Max)

*==>  Counters
      Integer  l,k
      Integer  I,J

      Integer Module, Tower, CAL, cType
      Integer get_fbr, get_mod, get_tow, get_typ
      Integer Cell ,Cal_CellNr
      Real    Cal_E,Cal_CellIm

*==>  define some usefull functions (F=0,B=1,R=2)
      get_fbr(Cell) = ISHFT(Cell,-14)
      get_Mod(Cell) = ISHFT(IAND(Cell,x'3E00'),-9)
      get_Tow(Cell) = ISHFT(IAND(Cell,x'1F0' ),-4)
      get_typ(Cell) = ISHFT(IAND(Cell,x'E'   ),-1)

*==>  Ini......................................................................
      Ierr   = 0
      Nzufos = 0

      Do I=1,nT_Max
       mTrI(I) = 0
       Do J=1,10
        wTrI(J,I) = 0
       EndDo
      EndDo

      Do I=1,nI_Max
       ITRm(I) = 0
       Do J=1,15
        ITrw(J,I) = 0
       EndDo
      EndDo

*==>  prepare info for CTD+CAL matching..........................................
      hnI = nIsl
      Do I=1,hnI
       emcI(I) = emcEIsl(I)
       heI( I) = eIsl(   I)
       hxI( I) = xIsl(   I)
       hyI( I) = yIsl(   I)
       hzI( I) = zIsl(   I)
       hrI( I) = rIsl(   I)
       hcI( I) = NrcIsl( I)
       hthI(I) = Atan2(SQRT(hxI(I)**2+hyI(I)**2),hzI(I)-Zvtx(3))
       hphI(I) = Atan2(hyI(I),hxI(I))
       If (hphI(I).LT.0.) hphI(I) = hphI(I) + TwoPi
      EndDo

*..................................... 
*-->   match coneIslands to Tracks
*.....................................
      Do l=1,hnI
       Coord(1)  = hxI( l)
       Coord(2)  = hyI( l)
       Coord(3)  = hzI( l)
       thIsl     = hthI(l)
       hEIsl     = heI( l)
       epsIsl    = emcI(l)/hEIsl
       DCoord(1) = 5.0
       DCoord(2) = 5.0
       DCoord(3) = 5.0

       DCA_hCut = Max(DCA_hCut_min,hrI(l))

       Do j=1,nT

*-->    Distance between track and cluster
        Call VCTDCA(vcthID(j),Coord,DCoord,DistCA,DCATerr,DCACerr)
        If (DistCA.GE.0.AND.DistCA.LE.DCA_hCut) Then

*-->     Island to Track match info
         If (ITrm(l).LT.15) Then
          ITrm(        l) = ITrm(l) + 1
          ITrw(ITrm(l),l) = j
         EndIf

*-->     Track to Island match info
         If (mTrI(j).LT.10) Then   
          mTrI(        j) = mTrI(j)+1
          wTrI(mTrI(j),j) = l 
         EndIf

        EndIf
       EndDo
      EndDo

*-->  pick up all the good Unmatched Tracks.........................
      Do j=1,nT
       If ( mTrI(j).EQ.0 .AND. qltTr(j).EQ.1 ) Then
        Tr_E = Sqrt(pTr(j)**2 + xpi2)
        Nzufos         = Nzufos + 1
        tufo(1,Nzufos) = 0
        tufo(2,Nzufos) = j
        tufo(3,Nzufos) = 0
        tufo(4,Nzufos) = 0
        zufo(1,Nzufos) = pTr(j)*Sin(thTr(j))*Cos(phTr(j))
        zufo(2,Nzufos) = pTr(j)*Sin(thTr(j))*Sin(phTr(j))
        zufo(3,Nzufos) = pTr(j)*Cos(thTr(j))
        zufo(4,Nzufos) = Tr_E        
       EndIf
      EndDo

*-->  Now mask out CAL energy due to tracking.......................
      Do l=1,hnI        
       E_Isl = heI(l)

       If (ITrm(l).GE.1) Then

        Do k=1,ITrm(l)        
         j     = ITrw(k,l)

         If (mTrI(j).EQ.1.AND.qltTr(j).EQ.1) Then
          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 10
          tufo(2,Nzufos) = j
          tufo(3,Nzufos) = l
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = pTr(j)*Sin(thTr(j))*Cos(phTr(j))
          zufo(2,Nzufos) = pTr(j)*Sin(thTr(j))*Sin(phTr(j))
          zufo(3,Nzufos) = pTr(j)*Cos(thTr(j))
          zufo(4,Nzufos) = Sqrt(pTr(j)**2 + xpi2)
      
          E_Isl = E_Isl - pTr(j)
         EndIf

        EndDo 
 
*==>    Pickup this Island if it is above a E_Isl_ThresHold...............
        If (E_Isl.GT.0.1) Then
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 32
         tufo(2,Nzufos) = 0
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = E_Isl*Sin(hthI(l))*Cos(hphI(l))
         zufo(2,Nzufos) = E_Isl*Sin(hthI(l))*Sin(hphI(l))
         zufo(3,Nzufos) = E_Isl*Cos(hthI(l)) 
         zufo(4,Nzufos) = E_Isl
        EndIf       

       Else                     !*-->  Unmatched Islands
        
        Nzufos         = Nzufos + 1
        tufo(1,Nzufos) = 31
        tufo(2,Nzufos) = 0
        tufo(3,Nzufos) = l
        tufo(4,Nzufos) = 0
        zufo(1,Nzufos) = E_Isl*Sin(hthI(l))*Cos(hphI(l))
        zufo(2,Nzufos) = E_Isl*Sin(hthI(l))*Sin(hphI(l))
        zufo(3,Nzufos) = E_Isl*Cos(hthI(l)) 
        zufo(4,Nzufos) = E_Isl
*     
       EndIf
*
      EndDo
*
      End
