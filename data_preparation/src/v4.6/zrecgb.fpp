      Subroutine z_RecGB(isdoTrk,isdoIsl,vtx,Ierr)
*
*===> USER HAS TO FILL Common /zIsl01/
*===> USER HAS TO FILL Common /zRec02/
*
*      INPUT:  isdoTrk  <--> =  1 CAL_or_CTD  zufos
*                            =  2 CAL_and_CTD zufos
*                            = -3 zufos=coneIslands
*              isdoIsl  <--> =0 do NOT build coneIslands
*                            =1 do     build coneIslands
*               vtx(3)  <--> Reconstructed vertex(1=x,2=y,3=z)
*
*     OUTPUT:
*
*       Ierr:   0  <--> Everything is cosher 
*              -1  <--> nTrks in VCTPAR > nT_Max
*
*  Description:
*
*   Common/zRec01/ 
*     Nzufos            <--> # of ZUFOs
*     tufo(4,0:Nzufos)  <--> Type of ZUFO information
*     zufo(4,0:Nzufos)  <--> 4 momentum of ZUFO
*                            (1=px,2=py,3=pz,4=E)
*
*   Common/zRec02/  
*     eNdis             <--> # of cells belonging to electron
*     ePrdis            <--> Probability of found electron
*     eEdis             <--> CAL Energy of electron
*     eXdis             <--> CAL X-position of electron
*     eYdis             <--> CAL Y-position of electron
*     eZdis             <--> CAL Z-position of electron
*
*   Common/zRec03/  
*     eTrMat            <--> # of Tracks which are matched 
*                            to electron within eCut_DCA 
*     eTrNear           <--> VCTPAR_ID of the nearest track
*     vcteID(0:eTrMat)  <--> VCTPAR_IDs of all the matched tracks to electron
*     eTrDCA            <--> Distance of Closest Approach 
*                            to electron by nearest track(eTrNear=VCTPAR_ID)
*   Common/zRec04/  
*     nT                <--> # of selected tracks for ZUFOs consideration
*                            (excluding electron tracks /zRec03/)
*     vcthID(0:nT)      <--> VCTPAR_ID of the track
*     qTr(   0:nT)      <--> charge of the track 
*     swmTr( 0:nT)      <--> SwimCode of the track
*     qltTr( 0:nT)      <--> quality of the track
*                            (=0 "bad"  trk DO NOT USE for ZUFOs,
*                             =1 "good" trk DO     USE for ZUFOs)
*                            quality is controled by sl_CUT,ptMin_CUT,ptMax_CUT
*     mTrI(  0:nT)      <--> # of Islands track has matched
*     pTr(   0:nT)      <--> Track momentum
*     dpTr(  0:nT)      <--> Error on track momentum
*     thTr(  0:nT)      <--> Theta of the track
*     phTr(  0:nT)      <--> Phi of the track
*     wTrI(0:mTrI,0:nT) <--> Index of the Island in Common/zIsl02/
*                            to which track was matched.
*
*   Author:   Gennady Briskin, Tel Aviv University
*   Date:     20-Dec-1996
*
*_______________________________________________________________________________
      Implicit NONE
*
#include "zisles.inc"
#include "zdrecgb.inc"

*-->  Input
      Integer isdoTrk,isdoIsl
      Real    vtx(3)

*==>  OutPut
      Integer Ierr

*==>  Counters
      Integer l

*==>  Useful Constants
      REAL       Pi, TwoPi
      PARAMETER (Pi = 3.14159265358979323846)
      PARAMETER (TwoPi = 2.*Pi)

*==>  Local variables
      Real    Zvtx(3)
      Real    thIsl,phIsl
      Real    xi, yi, zi

*==>  Ini things
      Ierr    = 0
      Nzufos  = 0

      Zvtx(1) = vtx(1)
      Zvtx(2) = vtx(2)
      Zvtx(3) = vtx(3)

*==>  Assume this event has no Tracks (zufos=coneIslands)
      If (isdoTrk.EQ.-3) nT = 0
      
*==>  get Tracking info...................................
      If (isdoTrk.GT.0) Then
       Call z_Trks(Ierr)
       If(Ierr.LT.0) Return
      EndIf

*==>  get coneIsland info.................................
      If (isdoIsl.EQ.1) Then
       Call Isles(Zvtx,PrCutHad,PrCutElec,Ierr)
      EndIf

*==>  And now lets have some fun %~{)
      If (nT.EQ.0) Then

*==>   take info from coneIslands.
       Nzufos = nIsl
       Do l=1,nIsl
        xi = xIsl(l)-Zvtx(1)
        yi = yIsl(l)-Zvtx(2)
        zi = zIsl(l)-Zvtx(3)
        
        thIsl = Atan2(Sqrt(xi**2+yi**2),zi)
        phIsl = Atan2(yi,xi)
        If (phIsl.LT.0.) phIsl = phIsl + TwoPi

        tufo(1,l) = 31
        tufo(2,l) = 0
        tufo(3,l) = l
        tufo(4,l) = 0
        zufo(1,l) = eIsl(l)*Sin(thIsl)*Cos(phIsl)
        zufo(2,l) = eIsl(l)*Sin(thIsl)*Sin(phIsl)
        zufo(3,l) = eIsl(l)*Cos(thIsl) 
        zufo(4,l) = eIsl(l)

       EndDo
       
      ElseIf(ABS(isdoTrk).EQ.1) Then

       Call z_CAL_or_CTD(Zvtx,Ierr)
 
      ElseIf(ABS(isdoTrk).EQ.2) Then

       Call z_CAL_and_CTD(Zvtx,Ierr)

      Else

       Write(6,*)'*==> z_RecGB: No such isdoTrk option ???',isdoTrk
       STOP       
       
      EndIf
*
*==>  Here goes the timing routine by AQ.

*
 123  FORMAT(2X,I5,1X,I6,1X,A50,I3)
      End

