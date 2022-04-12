      Subroutine z_CAL_or_CTD(Zvtx,Ierr)
      Implicit NONE
*----------------------------------------------------------------------
*
*   Called by z_RecGB routine to build CAL_or_CTD ZUFOs. 
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
*   --------------------------------------------------------------------
*   Author:  Niels Tuning -- February 1999 
*   Version: ZufosNT v4.3
*
*   Implemented changes made by Arnulf Quadt, Ulrike Wollmer and myself.
*   Documentation: www-zeus.desy.de/~ntuning/ZEUS_ONLY/zufo_new.html
*______________________________________________________________________
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
      Real    thIsl,hEIsl
      Real    matPr

      Real    Tr_E,E_Obj,E_Isl,epsIsl,P_Trs,dP_Trs
      Real    terr1,terr2
      Real    Cerr, DEovE, Dpovp, dEpx
      Real    DCA_hCut
      Real    conf_fac
      Logical L_take1to1_match,L_take2to1_match,L_take3to1_match,
     &        L_take1to2_match,L_take2to2_match,L_take1to2_trkp

      Integer j1,j2,j3,i1,i2,i3
      Real    t1,t2,t3
      Integer il11,il21,il12,il22

      Real    x1,x2,x3,y1,y2,y3,z1,z2,z3
      Real    ugl,xpipi2,xpipipi2

*==>  Island Info....................................................
      Integer  hnI         , hcI(   nI_Max)
      Real     heI( nI_Max), emcI(  nI_Max)
      Real     hxI( nI_Max), hyI(   nI_Max)
      Real     hzI( nI_Max), hrI(   nI_Max)
      Real     hthI(nI_Max), hphI(  nI_Max)
      Integer  ITRm(nI_Max),ITrw(15,nI_Max)

*==>  Counters........................................................
      Integer  l,k
      Integer  I,J

      Integer Module, Tower, CAL, cType
      Integer get_fbr, get_mod, get_tow, get_typ
      Integer Cell ,Cal_CellNr
      Real    Cal_E,Cal_CellIm

c --- secondary vertex tracks
      integer thisid


*==>  define some usefull functions (F=0,B=1,R=2).....................
      get_fbr(Cell) = ISHFT(Cell,-14)
      get_Mod(Cell) = ISHFT(IAND(Cell,x'3E00'),-9)
      get_Tow(Cell) = ISHFT(IAND(Cell,x'1F0' ),-4)
      get_typ(Cell) = ISHFT(IAND(Cell,x'E'   ),-1)

*==>  Ini..............................................................
      Ierr   = 0
      Nzufos = 0

      Do I=1,nT_Max
       mTrI(I)       = 0
       Dist_IsTr( I) = 999. 
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

*==>  prepare info for CTD+CAL matching.................................
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

*==>    Distance of Clossest Approach between track and cluster
          if (vcthID(j).gt.0) then
             Call VCTDCA(vcthID(j),Coord,DCoord,DistCA,DCATerr,DCACerr)
          elseif (vcthID(j).lt.0) then
             thisid = -vcthID(j)
             Call HLTDCA(thisid,Coord,DCoord,DistCA,DCATerr,DCACerr)
          endif
c         Call Matsch(vcthID(j),1,Coord,thIsl,hEIsl,epsIsl,1,vmode,
c     &               DistCA,DCACerr,matPr,Ierr)         
        If (DistCA.GE.0.AND.DistCA.LE.DCA_hCut) Then

*==>    NT; add DCA information to closest island for each trk 
           Dist_IsTr( j) = MIN(Dist_IsTr( j),DistCA)

*==>     Island to Track match info
         If (ITrm(l).LT.15) Then
          ITrm(        l) = ITrm(l) + 1
          ITrw(ITrm(l),l) = j
         EndIf

*==>     Track to Island match info
         If (mTrI(j).LT.10) Then   
          mTrI(        j) = mTrI(j)+1
          wTrI(mTrI(j),j) = l 
         EndIf

        EndIf
       EndDo
      EndDo

*==>  Now we deal with Unmatched Tracks.
      If (zRecGB_Debug) Then
       Write(6,603)
      EndIf
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
        
        If (zRecGB_Debug) Then
         Write(6,604)pTr(j),thTr(j)*180./Pi
        EndIf
       EndIf
      EndDo

*...........................................................
*-->  Now we calculate Next Generation hadronic variables
*...........................................................
      If (zRecGB_Debug) Then
       Write(6,600)
      EndIf

      Do l=1,hnI

*-->  Debug Print
       If (zRecGB_Debug) Then
        Write(6,601)l,heI(l),ITrm(l),hcI(l)
        Write(6,*)'----------------------------------------------'
        Write(6,605)
        Do k=1,hcI(l)
         Cal_CellNr = pCellIsl(  k,l) 
         Cal_E      = ECellIsl(  k,l)
         Cal_CellIm = ImbCellIsl(k,l)/Cal_E
         CAL        = get_fbr(Cal_CellNr)
         Module     = get_Mod(Cal_CellNr)+1
         Tower      = get_Tow(Cal_CellNr)+1          
         cType      = get_typ(Cal_CellNr)
         Write(6,606)k,CAL,Module,Tower,cType,Cal_E,Cal_CellIm
        EndDo
        Write(6,*)'----------------------------------------------'
        Do k=1,nT
         Do j=1,mTrI(k)
          If (wTrI(j,k).EQ.l) Then
           Write(6,602)k,l,mTrI(k),qltTr(k)
          EndIf
         EndDo
        EndDo
       EndIf
        
       If (ITrm(l).GE.1) Then
        
        j    = ITrw(1,l)
        Tr_E = Sqrt(pTr(j)**2 + xpi2)
        j1   = 0
        j2   = 0
        i1   = 0
        i2   = 0
        i3   = 0
        t1   = 0
        t2   = 0
        L_take1to1_match = .FALSE.
        L_take2to1_match = .FALSE.
        L_take3to1_match = .FALSE.
        L_take2to2_match = .FALSE.
        L_take1to2_trkp  = .FALSE.
c N.T. 2/99
        L_take1to2_match = .FALSE.
        

*-->    Don't assign mass before making some desicion first.
        E_Isl = -1.

*-->    matching logic
        If(ITrm(l).EQ.1.AND.mTrI(j).EQ.1)Then

         If (qltTr(j).EQ.1) Then
          E_Obj   = heI(l)   
          epsIsl = emcI(l)/E_Obj
          If    (epsIsl.GE.0.9) Then
           Cerr = 0.24
          Else
           Cerr = 0.4
          EndIf
          DEovE = Cerr/sqrt(E_Obj)
          Dpovp = dpTr(j)/pTr(j)
          terr1 = (Cerr*sqrt(E_Obj))/pTr(j)
          terr2 = (E_Obj*dpTr(j))/pTr(j)**2
          dEpx  = zNsig*Sqrt(terr1**2+terr2**2)
C ---     AQ+NT: in outer F/RCAL tracks favoured
          If ((hzI(l).lt.-140.0.OR.hzI(l).gt.220.0).AND.
     &        sqrt(hxI(l)**2+hyI(l)**2).gt.70.0    .AND.
     &        AQ_options) THEN    
            conf_fac=1.2
          ELSE
            conf_fac=1.0
          ENDIF

          If ((E_Obj/pTr(j)).LT.(EovPCut+dEpx).AND.
     &        Dpovp         .LT.DEovE*conf_fac) Then
           L_take1to1_match = .TRUE.
          ElseIf((E_Obj/pTr(j)).LT.(EovPCut+dEpx).AND.
     &           Dpovp         .GE.conf_fac*DEovE.AND.
     &           AQ_options) Then
           L_take1to2_trkp  = .TRUE.
          Else
*==>       This Island deserves at least a pion mass :-)
           E_Isl = Sqrt(E_Obj**2+xpi2)   
           If (zRecGB_Debug) Then
            Write(6,*)'*z_RecGB:1to1>',E_Obj,pTr(j),depx,DEovE,Dpovp
           EndIf
          EndIf

C ---    AQ: could this be a muon ?
         ElseIf (AQ_options.AND.
     &           qltTr(j).EQ.2) Then
          E_Obj   = heI(l)
          If ((E_Obj/pTr(j)).LT.muon_eovp.AND.
     &         E_Obj        .LT.muon_Emax) Then
           L_take1to1_match = .TRUE.
C ---     AQ: take island energy and track angles ...
          ElseIf((E_Obj/pTr(j)).LT.(EovPCut+dEpx).AND.
     &           Dpovp         .GE.conf_fac*DEovE.AND.
     &           AQ_options) Then
           L_take1to2_trkp  = .TRUE.
          Else
*==>       This Island deserves at least a pion mass :-)
           E_Isl = Sqrt(E_Obj**2+xpi2)
           If (zRecGB_Debug) Then
            Write(6,*)'*z_RecGB:1to1>',E_Obj,pTr(j),depx,DEovE,Dpovp
           EndIf
          EndIf
         EndIf

        ElseIf(ITrm(l).EQ.1.AND.mTrI(j).EQ.2)Then
         
         i1   = wTrI(1,j)
         i2   = wTrI(2,j)

c N.T. 2/99
        If( (qltTr(j).EQ.1.OR.qltTr(j).EQ.2).AND.
     &       ITrm(i1).EQ.1.AND.ITrm(i2).EQ.1) 
     &        L_take1to2_match = .TRUE.


        ElseIf(ITrm(l).EQ.1.AND.mTrI(j).EQ.3)Then

         i1   = wTrI(1,j)
         i2   = wTrI(2,j)
         i3   = wTrI(3,j)

        ElseIf (ITrm(l).EQ.2) Then

         j1   = ITrw(1,l)
         j2   = ITrw(2,l)
         t1   = sqrt(pTr(j1)**2 + xpi2)
         t2   = sqrt(pTr(j2)**2 + xpi2)
         
         If ( qltTr(j1).EQ.1.AND.qltTr(j2).EQ.1 )Then

*-->      both tracks matched this Island
          If ( mTrI(j1).EQ.1.AND. mTrI(j2).EQ.1 )Then
           E_Obj  = heI(l)   
           P_Trs  = pTr(j1)+pTr(j2)
           dP_Trs = sqrt(dpTr(j1)**2+dpTr(j2)**2)
           epsIsl = emcI(l)/E_Obj
           If    (epsIsl.GE.0.9) Then
            Cerr = 0.24
           Else
            Cerr = 0.4
           EndIf
           DEovE = Cerr/sqrt(E_Obj)
           Dpovp = dP_Trs/P_Trs
           terr1 = (Cerr*sqrt(E_Obj))/P_Trs
           terr2 = (E_Obj*dP_Trs)/P_Trs**2
           dEpx  = zNsig*Sqrt(terr1**2+terr2**2)
*     
           If ((E_Obj/P_Trs).LT.(EovPCut+dEpx).AND.Dpovp.LT.DEovE)Then
            L_take2to1_match = .TRUE.
           Else
*-->        This Island deserves at least a dipipi mass :-)
            ugl= Sin(thTr(j1))*sin(thTr(j2))*cos(phTr(j1)-phTr(j2))+
     +           cos(thTr(j1))*cos(thTr(j2))
            xpipi2=2.0*xpi2+2.0*(t1*t2 - pTr(j1)*pTr(j2)*ugl)
            E_Isl = Sqrt(E_Obj**2 + xpipi2)    
            
            If (zRecGB_Debug) Then
             Write(6,*)'*z_RecGB:2to1>',E_Obj,P_Trs,dEpx,DEovE,Dpovp
            EndIf
           EndIf

*-->      both tracks matched this island and some other ones          
          ElseIf (mTrI(j1).EQ.2.AND.mTrI(j2).EQ.2)Then
           il11 = wTrI(1,j1)
           il21 = wTrI(2,j1)
           il12 = wTrI(1,j2)
           il22 = wTrI(2,j2)
*-->       is this 2 to 2 match ?
           If(ITrm(il11).EQ.2.AND.ITrm(il21).EQ.2.AND.
     &          il11.EQ.il12.AND.il21.EQ.il22 )Then
*     
            E_Obj  = heI(il11)+heI(il22)   
            P_Trs  = pTr(j1)+pTr(j2)
            dP_Trs = sqrt(dpTr(j1)**2+dpTr(j2)**2)
            epsIsl = (emcI(il11)+emcI(il22))/E_Obj
            If    (epsIsl.GE.0.9) Then
             Cerr = 0.24
            Else
             Cerr = 0.4
            EndIf
            DEovE = Cerr/sqrt(E_Obj)
            Dpovp = dP_Trs/P_Trs
            terr1 = (Cerr*sqrt(E_Obj))/P_Trs
            terr2 = (E_Obj*dP_Trs)/P_Trs**2
            dEpx  = zNsig*Sqrt(terr1**2+terr2**2)
*     
            If ((E_Obj/P_Trs).LT.(EovPCut+dEpx).AND.Dpovp.LT.DEovE)Then
             L_take2to2_match = .TRUE.
            Else
             If (zRecGB_Debug) Then
              Write(6,*)'*z_RecGB:2to2>',E_Obj,P_Trs,depx,DEovE,Dpovp
             EndIf
            EndIf
            
           EndIf
          EndIf
         EndIf

        ElseIf (ITrm(l).EQ.3.OR.
     &          ITrm(l).GE.3.and.AQ_options) Then

         j1   = ITrw(1,l)
         j2   = ITrw(2,l)
         j3   = ITrw(3,l)
         t1   = sqrt(pTr(j1)**2 + xpi2)
         t2   = sqrt(pTr(j2)**2 + xpi2)
         t3   = sqrt(pTr(j3)**2 + xpi2)

         If(qltTr(j1).EQ.1.AND.qltTr(j2).EQ.1.AND.qltTr(j3).EQ.1)Then
          If(mTrI(j1).EQ.1.AND. mTrI(j2).EQ.1.AND. mTrI(j3).EQ.1)Then
           E_Obj  = heI(l)   
           P_Trs  = pTr(j1)+pTr(j2)+pTr(j3)
           dP_Trs = sqrt(dpTr(j1)**2+dpTr(j2)**2+dpTr(j3)**2)
           epsIsl = emcI(l)/E_Obj
           If    (epsIsl.GE.0.9) Then
            Cerr = 0.24
           Else
            Cerr = 0.4
           EndIf
           DEovE = Cerr/sqrt(E_Obj)
           Dpovp = dP_Trs/P_Trs
           terr1 = (Cerr*sqrt(E_Obj))/P_Trs
           terr2 = (E_Obj*dP_Trs)/P_Trs**2
           dEpx  = zNsig*Sqrt(terr1**2+terr2**2)
*     
           If ((E_Obj/P_Trs).LT.(EovPCut+dEpx).AND.Dpovp.LT.DEovE)Then
            L_take3to1_match = .TRUE.
           Else
*-->       This Island deserves at least a tripipipi mass :-)
            x1 = pTr(j1)*Sin(thTr(j1))*Cos(phTr(j1))
            x2 = pTr(j2)*Sin(thTr(j2))*Cos(phTr(j2))
            x3 = pTr(j3)*Sin(thTr(j3))*Cos(phTr(j3))
            y1 = pTr(j1)*Sin(thTr(j1))*Sin(phTr(j1))
            y2 = pTr(j2)*Sin(thTr(j2))*Sin(phTr(j2))
            y3 = pTr(j3)*Sin(thTr(j3))*Sin(phTr(j3))
            z1 = pTr(j1)*Cos(thTr(j1))
            z2 = pTr(j2)*Cos(thTr(j2))
            z3 = pTr(j3)*Cos(thTr(j3))
            xpipipi2 = (t1+t2+t3)**2
     -           - (x1+x2+x3)**2
     -           - (y1+y2+y3)**2
     -           - (z1+z2+z3)**2
            E_Isl    = Sqrt(E_Obj**2 + xpipipi2)   

            If (zRecGB_Debug) Then
             Write(6,*)'*z_RecGB:3to1>',E_Obj,P_Trs,depx,DEovE,Dpovp
            EndIf
           EndIf

          EndIf
         EndIf
        EndIf
*
        If (E_Isl.LT.0.) E_Isl = heI(l)

*-->    1 Track to 1 Island
        If(L_take1to1_match)Then
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 1
         tufo(2,Nzufos) = j
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = pTr(j)*Sin(thTr(j))*Cos(phTr(j))
         zufo(2,Nzufos) = pTr(j)*Sin(thTr(j))*Sin(phTr(j))
         zufo(3,Nzufos) = pTr(j)*Cos(thTr(j)) 
         zufo(4,Nzufos) = Tr_E
         
*-->    2 Tracks to 1 Island
        ElseIf(L_take2to1_match)Then
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 2
         tufo(2,Nzufos) = j1
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = pTr(j1)*Sin(thTr(j1))*Cos(phTr(j1))
         zufo(2,Nzufos) = pTr(j1)*Sin(thTr(j1))*Sin(phTr(j1))
         zufo(3,Nzufos) = pTr(j1)*Cos(thTr(j1)) 
         zufo(4,Nzufos) = t1
*
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 2
         tufo(2,Nzufos) = j2
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = pTr(j2)*Sin(thTr(j2))*Cos(phTr(j2))
         zufo(2,Nzufos) = pTr(j2)*Sin(thTr(j2))*Sin(phTr(j2))
         zufo(3,Nzufos) = pTr(j2)*Cos(thTr(j2)) 
         zufo(4,Nzufos) = t2

*-->    3 Tracks to 1 Island
        ElseIf(L_take3to1_match)Then
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 3
         tufo(2,Nzufos) = j1
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = pTr(j1)*Sin(thTr(j1))*Cos(phTr(j1))
         zufo(2,Nzufos) = pTr(j1)*Sin(thTr(j1))*Sin(phTr(j1))
         zufo(3,Nzufos) = pTr(j1)*Cos(thTr(j1)) 
         zufo(4,Nzufos) = t1
*     
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 3
         tufo(2,Nzufos) = j2
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = pTr(j2)*Sin(thTr(j2))*Cos(phTr(j2))
         zufo(2,Nzufos) = pTr(j2)*Sin(thTr(j2))*Sin(phTr(j2))
         zufo(3,Nzufos) = pTr(j2)*Cos(thTr(j2)) 
         zufo(4,Nzufos) = t2
*     
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 3
         tufo(2,Nzufos) = j3
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = pTr(j3)*Sin(thTr(j3))*Cos(phTr(j3))
         zufo(2,Nzufos) = pTr(j3)*Sin(thTr(j3))*Sin(phTr(j3))
         zufo(3,Nzufos) = pTr(j3)*Cos(thTr(j3)) 
         zufo(4,Nzufos) = t3
         
*-->    1 Track to 2 Islands(split Island case)
c        ElseIf(ITrm( l).EQ.1.AND.mTrI(j).EQ.2  .AND.
c     &         (qltTr(j).EQ.1.OR.qltTr(j).EQ.2).AND.
c     &        ITrm(i1).EQ.1.AND.ITrm(i2).EQ.1) Then
         ElseIf(L_take1to2_match) Then
c
c     N.T 2/99 (Rolf on mips3-sgi-irix6.2, 98a);
c Subscript out of range on file line 708, procedure z_cal_or_ctd.
c Attempt to access the 0-th element of variable ?.

C ---            this case is delt with later ...

  
*-->    1 Track to 3 Islands(split Island case)
c        ElseIf(ITrm( l).EQ.1.AND.mTrI(j).EQ.3.AND.qltTr(j).EQ.1.AND.
c     &         ITrm(i1).EQ.1.AND.
c     &         ITrm(i2).EQ.1.AND.
c     &         ITrm(i3).EQ.1) Then
c
c               WRITE(*,*) 'dubious case ... !!!'
C ---            this case is delt with later ...


*-->    2 Tracks to 2 Islands
        ElseIf(L_take2to2_match)Then
         If (l.EQ.il22) Then
          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 22
          tufo(2,Nzufos) = j1
          tufo(3,Nzufos) = il11
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = pTr(j1)*Sin(thTr(j1))*Cos(phTr(j1))
          zufo(2,Nzufos) = pTr(j1)*Sin(thTr(j1))*Sin(phTr(j1))
          zufo(3,Nzufos) = pTr(j1)*Cos(thTr(j1)) 
          zufo(4,Nzufos) = t1
*
          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 22
          tufo(2,Nzufos) = j2
          tufo(3,Nzufos) = il21
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = pTr(j2)*Sin(thTr(j2))*Cos(phTr(j2))
          zufo(2,Nzufos) = pTr(j2)*Sin(thTr(j2))*Sin(phTr(j2))
          zufo(3,Nzufos) = pTr(j2)*Cos(thTr(j2)) 
          zufo(4,Nzufos) = t2
*
          If (zRecGB_Debug) Then
           Write(6,*)'2to2>Isl',l,' matched by',ITrm(l),' Trs',j1,j2
          EndIf
         Else
          If (zRecGB_Debug) Then
           Write(6,*)'>Isl',l,heI(l),' matched by',ITrm(l),' Trs',j1,j2
           Write(6,*)'>Isl',il22,' matched by',ITrm(il22),' Trs',j1,j2
           Write(6,*)'>Trk',j1,pTr(j1),' Islands',il11,heI(il11),
     &          il21,heI(il21)
           Write(6,*)'>Trk',j2,pTr(j2),' Islands',il12,heI(il12),
     &          il22,heI(il22)
          EndIf
         EndIf            

C ---   AQ: track match, island taken for ebergy, track for position
        ElseIf(L_take1to2_trkp)Then
         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 41
         tufo(2,Nzufos) = j
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = heI(l)*Sin(thTr(j))*Cos(phTr(j))
         zufo(2,Nzufos) = heI(l)*Sin(thTr(j))*Sin(phTr(j))
         zufo(3,Nzufos) = heI(l)*Cos(thTr(j))
         zufo(4,Nzufos) = heI(l)


*-->     Some other configuration: track match, but island taken
        Else
         
C ---    AQ: Do we always want to use the island or 
C            do we want to use track sometimes ... ???

         Nzufos         = Nzufos + 1
         tufo(1,Nzufos) = 30
         tufo(2,Nzufos) = 0
         tufo(3,Nzufos) = l
         tufo(4,Nzufos) = 0
         zufo(1,Nzufos) = heI(l)*Sin(hthI(l))*Cos(hphI(l))
         zufo(2,Nzufos) = heI(l)*Sin(hthI(l))*Sin(hphI(l))
         zufo(3,Nzufos) = heI(l)*Cos(hthI(l)) 
         zufo(4,Nzufos) = E_Isl
         
        EndIf
        
       Else                     !*-->  Unmatched Islands
        
        Nzufos         = Nzufos + 1
        tufo(1,Nzufos) = 31
        tufo(2,Nzufos) = 0
        tufo(3,Nzufos) = l
        tufo(4,Nzufos) = 0
        zufo(1,Nzufos) = heI(l)*Sin(hthI(l))*Cos(hphI(l))
        zufo(2,Nzufos) = heI(l)*Sin(hthI(l))*Sin(hphI(l))
        zufo(3,Nzufos) = heI(l)*Cos(hthI(l)) 
        zufo(4,Nzufos) = heI(l)
*     
       EndIf
*
      EndDo

*-->  1 Track to 2 Islands
      Do j=1,nT
       If (mTrI(j).EQ.2.AND.(qltTr(j).EQ.1.OR.qltTr(j).EQ.2))Then
        i1   = wTrI(1,j)
        i2   = wTrI(2,j)        
        If(ITRm(i1).EQ.1.AND.ITRm(i2).EQ.1)Then
*     
         Tr_E   = Sqrt(pTr(j)**2 + xpi2)
         E_Obj  = heI(i1)+heI(i2)   
         epsIsl = (emcI(i1)+emcI(i2))/E_Obj
         If    (epsIsl.GE.0.9) Then
          Cerr = 0.24
         Else
          Cerr = 0.4
         EndIf
         DEovE = Cerr/sqrt(E_Obj)
         Dpovp = dpTr(j)/pTr(j)
         terr1 = (Cerr*sqrt(E_Obj))/pTr(j)
         terr2 = (E_Obj*dpTr(j))/pTr(j)**2
         dEpx  = zNsig*Sqrt(terr1**2+terr2**2)

C ---    AQ: in outer F/RCAL tracks favoured
         If ((hzI(i1).lt.-140.0.OR.hzI(i1).gt.220.0).AND.
     &       sqrt(hxI(i1)**2+hyI(i1)**2).gt.70.0    .AND.
     &       (hzI(i2).lt.-140.0.OR.hzI(i2).gt.220.0).AND.
     &       sqrt(hxI(i2)**2+hyI(i2)**2).gt.70.0    .AND.
     &       AQ_options) THEN    
           conf_fac=1.2
         ELSE
           conf_fac=1.0
         ENDIF
*
         L_take1to2_match  = .FALSE.
         L_take1to2_trkp   = .FALSE.
         If ((E_Obj/pTr(j)).LT.(EovPCut+dEpx) .AND.
     &       Dpovp         .LT.conf_fac*DEovE .AND.
     &       qltTr(j)      .EQ.1) Then
          L_take1to2_match = .TRUE.
         ElseIf((E_Obj/pTr(j)).LT.(EovPCut+dEpx)     .AND.
     &          Dpovp         .GE.conf_fac*DEovE     .AND.
     &          (qltTr(j)     .EQ.1.OR.qltTr(j).EQ.2).AND.
     &          AQ_options) Then
          L_take1to2_trkp  = .TRUE.
         Else
          If (zRecGB_Debug) Then
           Write(6,*)'*z_RecGB:1to2',E_Obj,pTr(j),depx,DEovE,Dpovp
          EndIf
         EndIf
*     
         If(L_take1to2_match)Then
          
          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 12
          tufo(2,Nzufos) = j
          tufo(3,Nzufos) = i1
          tufo(4,Nzufos) = i2
          zufo(1,Nzufos) = pTr(j)*Sin(thTr(j))*Cos(phTr(j))
          zufo(2,Nzufos) = pTr(j)*Sin(thTr(j))*Sin(phTr(j))
          zufo(3,Nzufos) = pTr(j)*Cos(thTr(j)) 
          zufo(4,Nzufos) = Tr_E
*     
C ---    AQ: take CAL energy and track position ---
         ElseIf(L_take1to2_trkp)Then

          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 37
          tufo(2,Nzufos) = j
          tufo(3,Nzufos) = i1
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = heI(i1)*Sin(thTr(j))*Cos(phTr(j))
          zufo(2,Nzufos) = heI(i1)*Sin(thTr(j))*Sin(phTr(j))
          zufo(3,Nzufos) = heI(i1)*Cos(thTr(j))
          zufo(4,Nzufos) = heI(i1)
*
          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 37
          tufo(2,Nzufos) = j
          tufo(3,Nzufos) = i2
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = heI(i2)*Sin(thTr(j))*Cos(phTr(j))
          zufo(2,Nzufos) = heI(i2)*Sin(thTr(j))*Sin(phTr(j))
          zufo(3,Nzufos) = heI(i2)*Cos(thTr(j))
          zufo(4,Nzufos) = heI(i2)
*
         Else

C --- AQ  Don't we want to match the track to one object
C         and use the cal only for the other one ... ???

          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 30
          tufo(2,Nzufos) = 0
          tufo(3,Nzufos) = i1
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = heI(i1)*Sin(hthI(i1))*Cos(hphI(i1))
          zufo(2,Nzufos) = heI(i1)*Sin(hthI(i1))*Sin(hphI(i1))
          zufo(3,Nzufos) = heI(i1)*Cos(hthI(i1))
          zufo(4,Nzufos) = heI(i1)
*
          Nzufos         = Nzufos + 1
          tufo(1,Nzufos) = 30
          tufo(2,Nzufos) = 0
          tufo(3,Nzufos) = i2
          tufo(4,Nzufos) = 0
          zufo(1,Nzufos) = heI(i2)*Sin(hthI(i2))*Cos(hphI(i2))
          zufo(2,Nzufos) = heI(i2)*Sin(hthI(i2))*Sin(hphI(i2))
          zufo(3,Nzufos) = heI(i2)*Cos(hthI(i2))
          zufo(4,Nzufos) = heI(i2)
*
         EndIf
        EndIf
       EndIf
      EndDo


*     
 600  Format('*z_RecGB>   Isl  Energy_Isl  matbyTrks   #Cells')
 601  Format(   12X,      I3,2X,F9.4   ,3X,I7,3X,I6)
 602  Format(3X,'*Trk->',I3,' --->',I2,' and has',I3,' matche(s)')
 603  Format('*z_RecGB> NOT matched TrP  Theta')
 604  Format(10X,F9.3,2X,F5.1)
 605  Format(' Cell  CAL  Module  Tower  Type   Ecell  Imbalance')
 606  Format(1X,I4,2X,I3,2X,I6,2X,I5,2X,I4,2X,F7.3,1X,F8.2)
*
      End

