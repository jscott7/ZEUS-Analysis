      Subroutine Isles(Vertex,ProbCuth,ProbCute,Ierr)
*
*
* This routine performs CAL clustering by matching cell Islands
* based on angular separation. The Cell Islands are found for each
* CAL separately (F,R,B), and separately for each section (EMC,HAC1,HAC2)
* using the cell Islands routine developed by Gennady Briskin.  (The original
* Island algorithm was developed by Larry Wai.)  A probability is assigned
* based on the 3-D angular separation from on a parametrization derived
* from a single-pion MC.
*
* USER INPUT:  Vertex(3) -- Reconstructed vertex(1=x,2=y,3=z)
*
*              ProbCuth  -- probability cut for matching HAC2 -> HAC1
*                                                    and HAC1 -> EMC
*                           recommended value is 0.01 for general purpose,
*                           but can be looser in clean events (e.g., VM)
*
*              ProbCute  -- probability cut for matching EMC -> EMC
*                           recommended value is 1.0 (no matching allowed)
*                           for electrons, 0.02 for hadrons
*                           but can be looser in clean events (e.g., VM)
*
* OUTPUT:      Ierr    0 -- Everything is fine
*                      1 -- No CAL information (zCells=0)
*
* +CDE, ZIsles:
*
*  Isles Input:
*
*===>  COMMON /zIsl01/ HAS TO BE FILLED BY USER BEFORE CALLING ISLES ROUTINE 
*
*      COMMON /zIsl01/ 
*      zCells              -- Number of cells to be considered for clustering.
*      zPnrl(   1:Nr_cMax) -- PoserNr   of the cell
*      zEl(     1:Nr_cMax) -- Energy    of the cell
*      zImbl(   1:Nr_cMax) -- Imbalance of the cell
*      zT1(     1:Nr_cMax) -- time(1)   of the cell
*      zT2(     1:Nr_cMax) -- time(2)   of the cell
*      zID(     1:Nr_cMax) -- Row number in table Caltru(i.e. Caltru_ID)
*
*  Isles Output:
*
*      COMMON /zIsl02/ 
*
*      COMMON /zIsl03/
* 
*      COMMON /zIsl04/ 
*      I_caltru(1:Nr_cMax) -- cluster assignment for each cell in 
*                             COMMON /zIsl01/
*
*
* A. Caldwell 22/4/96
*
* -- modified --
* G. Briskin 24/10/97
*
******************************************************************************
*
      Implicit NONE
*
#include "zisles.inc"
!  /*       <-- Islands information  */
*
*==>  Input Variables
      Real    Vertex(3)
      Real    ProbCuth,ProbCute

*==>  Output
      Integer Ierr

*------------------------------
*     Internal Variables      *
*------------------------------

*==>  Counter for coneIslands
      Integer nclus

*==>  Some Usefull Constants
      REAL       Pi, TwoPi, PiBy2
      PARAMETER (Pi = 3.14159265358979323846)
      PARAMETER (TwoPi = 2.*Pi)
      PARAMETER (PiBy2 = Pi/2.0) 
      Real       todeg
      Parameter (todeg = 180./Pi)

*==>  Error flags
      Logical Lerr

*==>  Counters
      Integer Icell,no_I,n,m,l,I,K
      Integer kIsl,NrCells

*==>  Position recon variables
      Real    xC,yC,zC
      Real    dx, Er, El, d1, d2(2), wi, wtot

*==>  Isles variables
      Integer NEMC,NHAC1,NHAC2

      Real    E_EMC(   nI_Max),theta(    nI_Max),phi(      nI_Max)
      Real    E_HAC1(  nI_Max),theta1(   nI_Max),phi1(     nI_Max)
      Real    E_HAC2(  nI_Max),theta2(   nI_Max),phi2(     nI_Max)

      Integer EMCtoIsl(nI_Max),H1toIsl(  nI_Max),H2toIsl(  nI_Max)
      Integer match(   nI_Max),match1(   nI_Max),match2(   nI_Max)
      Integer nassign( nI_Max),nassignh1(nI_Max),nassignh2(nI_Max)
      Integer IND(     nI_Max)

*==>  Auxilary variables for angle determination
      Real    r,z_e

      Real    probmin,dphi12,dcos,dist,prob
      Real    probisles

*==>  Useful variables
      Integer CellNr,i_CAL,i_MOD,i_TOW,i_Typ

*==>  Variables for projection on perpendicular plane
      Real    rayX,rayY,rayZ
      Real    lenC,unit_Cx,unit_Cy,unit_Cz
      Real    rhoX,rhoY,rhoZ
      Real    Swim

*==>  Function Declarations
      Integer get_fbr, get_mod, get_tow, get_typ

      get_fbr(CellNr) = ISHFT(CellNr,-14)+1
      get_Mod(CellNr) = ISHFT(IAND(CellNr,x'3E00'),-9)+1
      get_Tow(CellNr) = ISHFT(IAND(CellNr,x'1F0' ),-4)+1
      get_typ(CellNr) = ISHFT(IAND(CellNr,x'E'   ),-1)

*==>  Zero Ground
      Ierr  = 0
      nIsl  = 0
      Nclus = 0

*==>  Check that we dont have too many or too few cells
      If (zCells.LE.0) then
       Ierr = 1
       Return
      Endif

*==> Use cellIslands routine from Gennady to get Islands
*    (mode 11=corner cells connect ,12=corner cells DONT connect)
      
      NEMC  = 0
      NHAC1 = 0 
      NHAC2 = 0 

      call z_cIsles(Vertex,zCells,zEl,zPnrl,zImbl,no_I,zcIsl_Mode,Ierr)
      If (Ierr.LT.0) then
       Write(6,*)'*==> Isles: Error from z_cIsles',Ierr
       Ierr = 2
       Return
      Endif

*==>  Get theta and phi of cellIslands
      Do kIsl=1,no_I

       If    (IslTyp(kIsl).EQ.1000) Then

        NEMC           = NEMC + 1

*==>    Store Island number +1000 for Caltru row number
        NrCells = NrcIsl(kIsl)
        Do I=1,NrCells
         l = zIsl_Cell(I,kIsl)
         I_Caltru(l)=NEMC+1000
        EndDo

        EMCtoIsl(NEMC) = kIsl
        E_EMC( NEMC)   = eIsl(kIsl)
        z_e            = zIsl(kIsl)-Vertex(3)
        r              = SQRT(yIsl(kIsl)**2+xIsl(kIsl)**2)
        theta( NEMC)   = ATAN2(r,z_e)
        phi(   NEMC)   = ATAN2(yIsl(kIsl),xIsl(kIsl))
        If(phi(NEMC).LT.0) phi(NEMC) = phi(NEMC) + TwoPi

        If (Isles_Debug) Then 
         Print *,'*--> EMC Island ',kIsl,eIsl(kIsl),NrcIsl(kIsl),
     &        xIsl(kIsl),yIsl(kIsl),zIsl(kIsl)
         NrCells = NrcIsl(kIsl)
         Do I=1,NrCells
          CellNr = pCellIsl(I,kIsl)
          i_CAL = get_fbr(CellNr)
          i_MOD = get_mod(CellNr)
          i_TOW = get_tow(CellNr)
          i_Typ = get_typ(CellNr)
          Print *,'*-> Cell,CAL,MOD,TOW,Type',
     &         CellNr,i_CAL,i_MOD,i_TOW,i_Typ
         EndDo
        EndIf

       ElseIf(IslTyp(kIsl).EQ.2000) Then

        NHAC1          = NHAC1 + 1

*==>    Store Island number +2000 for Caltru row number
        NrCells = NrcIsl(kIsl)
        Do I=1,NrCells
         l = zIsl_Cell(I,kIsl)
         I_Caltru(l)=NHAC1+2000
        EndDo

        H1toIsl(NHAC1) = kIsl
        E_HAC1( NHAC1) = eIsl(kIsl)
        z_e            = zIsl(kIsl)-Vertex(3)
        r              = SQRT(yIsl(kIsl)**2+xIsl(kIsl)**2)
        theta1( NHAC1) = ATAN2(r,z_e)
        phi1(   NHAC1) = ATAN2(yIsl(kIsl),xIsl(kIsl))
        If(phi1(NHAC1).LT.0) phi1(NHAC1) = phi1(NHAC1) + TwoPi

        If (Isles_Debug) Then
         Print *,'*--> HAC1 Island ',kIsl,eIsl(kIsl),NrcIsl(kIsl),
     &        xIsl(kIsl),yIsl(kIsl),zIsl(kIsl)
         NrCells = NrcIsl(kIsl)
         Do I=1,NrCells
          CellNr = pCellIsl(I,kIsl)
          i_CAL = get_fbr(CellNr)
          i_MOD = get_mod(CellNr)
          i_TOW = get_tow(CellNr)
          i_Typ = get_typ(CellNr)
          Print *,'*-> Cell,CAL,MOD,TOW,Type',I,i_CAL,i_MOD,i_TOW,i_Typ
         EndDo
        EndIf

       ElseIf(IslTyp(kIsl).EQ.3000) Then
 
        NHAC2          = NHAC2 + 1

*==>    Store Island number +3000 for Caltru row number
        NrCells = NrcIsl(kIsl)
        Do I=1,NrCells
         l = zIsl_Cell(I,kIsl)
         I_Caltru(l)=NHAC2+3000
        EndDo

        H2toIsl(NHAC2) = kIsl
        E_HAC2( NHAC2) = eIsl(kIsl)
        z_e            = zIsl(kIsl) - Vertex(3)
        r              = SQRT(yIsl(kIsl)**2+xIsl(kIsl)**2)
        theta2( NHAC2) = ATAN2(r,z_e)
        phi2(   NHAC2) = ATAN2(yIsl(kIsl),xIsl(kIsl))
        If(phi2(NHAC2).LT.0) phi2(NHAC2) = phi2(NHAC2) + TwoPi

        If (Isles_Debug) Then
         Print *,'*--> HAC2 Island ',kIsl,eIsl(kIsl),NrcIsl(kIsl),
     &   xIsl(kIsl),yIsl(kIsl),zIsl(kIsl)
         NrCells = NrcIsl(kIsl)
         Do I=1,NrCells
          CellNr = pCellIsl(I,kIsl)
          i_CAL = get_fbr(CellNr)
          i_MOD = get_mod(CellNr)
          i_TOW = get_tow(CellNr) 
          i_Typ = get_typ(CellNr)
          Print *,'*-> Cell,CAL,MOD,TOW,Type',I,i_CAL,i_MOD,i_TOW,i_Typ
         EndDo
        EndIf

       Else
        Print *,'*-> Unknown Island Type ???',kIsl,IslTyp(kIsl)
        STOP
       EndIf

      Enddo

************************************************************************
*  Lets do the matching
*
*  We start from hac2 and work our way inwards to the EMC.  Every EMC
* Island will be a seed for a 'super'-Island.  I.e., there is initially
* no grouping of EMC Islands.  EMCs can be matched if they are close 
* together depending on user supplied probability cut.
*
*  Get distances between clusters in theta-phi space and take the highest 
* probability match.  It has to be larger than the user supplied cut.
*************************************************************************

*
*==>  Start with HAC2
      Do n=1,nhac2
       match2(n) = 0
       probmin   = probcuth
*
*==>   Check for HAC1 within angular cone
       Do m=1,nhac1
        dphi12 = phi2(n)-phi1(m)
        dcos   = sin(theta2(n))*sin(theta1(m))*cos(dphi12)+
     &           cos(theta2(n))*cos(theta1(m))
        dist   = acos(dcos)

*==>    Get probability for this angular separation
        prob   = probisles(dist,2)

*==>    Keep highest probability match above cut
        If (prob.gt.probmin) then
         probmin   = prob
         match2(n) = m

         If (Isles_Debug) Then
          Print *,'*'
          Print *,'*'
          Print *,'*-> HAC2',n,' ->HAC1',m,' Prob',prob
          kIsl = H2toIsl(n)
          NrCells = NrcIsl(kIsl)
          Print *,'*'
          Print *,'*-> HAC2 Island: ',theta2(n)*todeg,phi2(n)*todeg
          Do I=1,NrCells
           CellNr = pCellIsl(I,kIsl)
           i_CAL = get_fbr(CellNr)
           i_MOD = get_mod(CellNr)
           i_TOW = get_tow(CellNr) 
           i_Typ = get_typ(CellNr)
           Print *,'*-> Cell,CAL,MOD,TOW,Type',I,i_CAL,i_MOD,i_TOW,i_Typ
          EndDo
*
          kIsl = H1toIsl(m)
          NrCells = NrcIsl(kIsl)
          Print *,'*'
          Print *,'*-> HAC1 Island',theta2(m)*todeg,phi1(m)*todeg
          Do I=1,NrCells
           CellNr = pCellIsl(I,kIsl)
           i_CAL = get_fbr(CellNr)
           i_MOD = get_mod(CellNr)
           i_TOW = get_tow(CellNr) 
           i_Typ = get_typ(CellNr)
           Print *,'*-> Cell,CAL,MOD,TOW,Type',I,i_CAL,i_MOD,i_TOW,i_Typ
          EndDo
         EndIf

        Endif
       EndDo
*
*==>  If there's no link to a HAC1, try a direct link to an EMC
       Do m=1,nemc
        dphi12 = phi2(n)-phi(m)
        dcos   = sin(theta2(n))*sin(theta(m))*cos(dphi12)+
     &           cos(theta2(n))*cos(theta(m))
        dist   = acos(dcos)
        prob   = probisles(dist,1)

*==>    Keep highest probability match above cut
        If (prob.gt.probmin) then
         probmin=prob
         match2(n) = m+1000

         If (Isles_Debug)
     &        print *,'hac2 -> EMC link',n,' emc',m,' prob',prob

        Endif
       EndDo
      EndDo
*
*==>  Now link the hac1 to emc
      Do n=1,nhac1
       probmin=probcuth
       match1(n)=0
       Do m=1,nemc
        dphi12 = phi1(n)-phi(m)
        dcos=sin(theta1(n))*sin(theta(m))*cos(dphi12)+
     &       cos(theta1(n))*cos(theta(m))
        dist=acos(dcos)
        prob = probisles(dist,1)
        If (prob.gt.probmin) then
         probmin=prob
         match1(n) = m

         If (Isles_Debug) Then
          Print *,'*'
          Print *,'*'
          Print *,'*-> HAC1',n,' ->EMC',m,' Prob',prob
          kIsl = H1toIsl(n)
          NrCells = NrcIsl(kIsl)
          Print *,'*'
          Print *,'*-> HAC1 Island: ',theta1(n)*todeg,phi1(n)*todeg
          Do I=1,NrCells
           CellNr = pCellIsl(I,kIsl)
           i_CAL = get_fbr(CellNr)
           i_MOD = get_mod(CellNr)
           i_TOW = get_tow(CellNr) 
           i_Typ = get_typ(CellNr)
           Print *,'*-> Cell,CAL,MOD,TOW,Type',I,i_CAL,i_MOD,i_TOW,i_Typ
          EndDo
*
          kIsl = EMCtoIsl(m)
          NrCells = NrcIsl(kIsl)
          Print *,'*'
          Print *,'*-> EMC Island',theta(m)*todeg,phi(m)*todeg
          Do I=1,NrCells
           CellNr = pCellIsl(I,kIsl)
           i_CAL = get_fbr(CellNr)
           i_MOD = get_mod(CellNr)
           i_TOW = get_tow(CellNr) 
           i_Typ = get_typ(CellNr)
           Print *,'*-> Cell,CAL,MOD,TOW,Type',I,i_CAL,i_MOD,i_TOW,i_Typ
          EndDo
         EndIf

        Endif
       enddo
      enddo
*
*  Now link the emc to emc
*
*  First sort in energy (ascending energy)
*
      CALL SORTZV(E_EMC,IND,nemc,1,0,0)
*     
      do n=1,nemc
*
*  Different probability requirement for matching EMC's
*
       probmin=probcute
       match(IND(n))=0
       do m=n+1,nemc
        dphi12 = phi(IND(n))-phi(IND(m))
        dcos=sin(theta(IND(n)))*sin(theta(IND(m)))*cos(dphi12)+
     &       cos(theta(IND(n)))*cos(theta(IND(m)))
        dist=acos(dcos)
        prob = probisles(dist,1)
        If (prob.gt.probmin) then
         probmin=prob
         match(IND(n)) = IND(m)

         If (Isles_Debug) Then
          Print *,'*' 
          Print *,'*' 
          Print *,'*-> EMC',IND(n),' -->EMC',IND(m),' prob',prob
          kIsl = EMCtoIsl(n)
          NrCells = NrcIsl(kIsl)
          Print *,'*'
          Print *,'*-> EMC Island: ',IND(n),theta(n)*todeg,phi(n)*todeg
          Do I=1,NrCells
           CellNr = pCellIsl(I,kIsl)
           i_CAL = get_fbr(CellNr)
           i_MOD = get_mod(CellNr)
           i_TOW = get_tow(CellNr) 
           i_Typ = get_typ(CellNr)
           Print *,'*-> Cell,CAL,MOD,TOW,Type',
     &          CellNr,i_CAL,i_MOD,i_TOW,i_Typ
          EndDo
*     
          kIsl = EMCtoIsl(m)
          NrCells = NrcIsl(kIsl)
          Print *,'*'
          Print *,'*-> EMC Island: ',IND(m),theta(m)*todeg,phi(m)*todeg
          Do I=1,NrCells
           CellNr = pCellIsl(I,kIsl)
           i_CAL = get_fbr(CellNr)
           i_MOD = get_mod(CellNr)
           i_TOW = get_tow(CellNr) 
           i_Typ = get_typ(CellNr)
           Print *,'*-> Cell,CAL,MOD,TOW,Type',
     &          CellNr,i_CAL,i_MOD,i_TOW,i_Typ
          EndDo
         EndIf

        Endif
       enddo
      enddo

*********************************************************************
* Now we update the cluster information
* We assign a cluster number to each CALTRU row.
*********************************************************************
*
* Zero the cluster assignment array for HAC2
*
      Do n=1,nhac2
       nassignh2(n) = 0
      Enddo
*
*  First assign the HAC2 cluster numbers to HAC1s
*
* Loop over caltru rows
      Do n=1,zCells
* Check is HAC2
       If (I_caltru(n).ge.3000) then
* If so, get cluster index
        m=I_caltru(n)-3000
* Check if this HAC2 island is matched to a HAC1 island
        If (match2(m).gt.0.and.match2(m).lt.1000) then
* If so, then point this caltru row to the HAC1 cluster
         I_caltru(n) = 2000+match2(m)
* Check if this HAC2 island is matched to an EMC island
        ElseIf (match2(m).gt.1000) then
* If so, then point this caltru row to the EMC cluster
         I_caltru(n) = match2(m)
        Else
* If not, its a standalone cluster
* Check if this cluster has a cluster number
         If (nassignh2(m).eq.0) then
* No, so give it one
          nclus = nclus + 1
          nassignh2(m) = nclus
         Endif
         I_caltru(n) = nassignh2(m)
        Endif
       Endif
      Enddo
*
*  Now HAC1s to EMCs
*
* Zero the cluster assignment array for HAC1
*
      Do n=1,nhac1
       nassignh1(n) = 0
      Enddo
*
      Do n=1,zCells
       If (I_caltru(n).ge.2000) then
        m=I_caltru(n)-2000
        If (match1(m).gt.0) then
         I_caltru(n) = 1000+match1(m)
        Else
* Check if this cluster has a cluster number
         If (nassignh1(m).eq.0) then
* No, so give it one
          nclus = nclus + 1
          nassignh1(m) = nclus
         Endif
         I_caltru(n) = nassignh1(m)
        Endif
       Endif
      Enddo
*
*  Now EMCs to EMCs -- note that their could be a chain n -> m -> k ...
*
      Do n=1,zCells
C	 print *,I_caltru(n)
       If (I_caltru(n).ge.1000) then
 19     m=I_caltru(n)-1000
C     Print *,m,match(m)
        If (match(m).gt.0) then
         I_caltru(n) = 1000+match(m)
        Else
         goto 20
        Endif
        goto 19
       Endif
 20    continue
      Enddo
*
* Zero the cluster assignment array
*
      Do n=1,nI_Max
       nassign(n) = 0
      Enddo
*
      Do n=1,zCells
       If (I_caltru(n).ge.1000) then
        m = I_caltru(n) - 1000
* Check if this cluster has a cluster number
        If (nassign(m).eq.0) then
*     No, so give it one
         nclus = nclus + 1
         nassign(m) = Nclus
        Endif
        I_caltru(n) = nassign(m)
       Endif
      Enddo

*==>  How many coneIslands did we find :-)
      nIsl = Nclus

**********************************************************************
*  We're done - now for some bookkeeping
**********************************************************************

*
* We now loop over the CALTRU cells and fill Gennady's arrays
*
      Do n=1,nI_Max
       eIsl(   n) = 0.
       xIsl(   n) = 0.
       yIsl(   n) = 0.
       zIsl(   n) = 0.
       rIsl(   n) = 0.
       xI1(    n) = 0.
       yI1(    n) = 0.
       zI1(    n) = 0.
       emcEIsl(n) = 0.
       NrcIsl( n) = 0
      Enddo
*     
      Do Icell=1,zCells
       kIsl = I_caltru(Icell)
       If (kIsl.GT.0.AND.kIsl.LE.nIsl) Then

        NrcIsl(                 kIsl) = NrcIsl(kIsl) + 1
*
        zIsl_Cell( NrcIsl(kIsl),kIsl) =       Icell
*
        ECellIsl(  NrcIsl(kIsl),kIsl) = zEl(  Icell) 
        ImbCellIsl(NrcIsl(kIsl),kIsl) = zImbl(Icell) 
        pCellIsl(  NrcIsl(kIsl),kIsl) = zPnrl(Icell) 
        eIsl(                   kIsl) = eIsl(kIsl) + zEl(Icell)
*        
        If (Isles_Debug) Then
         CellNr= zPnrl(Icell) 
         i_CAL = get_fbr(CellNr)
         i_MOD = get_mod(CellNr)
         i_TOW = get_tow(CellNr)
         i_Typ = get_typ(CellNr)
         Print *,'*-> Cell #',CellNr,'  Assigned to cluster ',kIsl
         Print *,'*-> CAL,MOD,TOW,Type',i_CAL,i_MOD,i_TOW,i_Typ
        Endif

       Else
        Write(6,*),'Isles: wrong # Isls ',kIsl,zCells,NEMC,NHAC1,NHAC2
        STOP
       Endif
      Enddo	 
*
*-->  get Island position
*
      Do kIsl=1,nIsl
*
       NrCells = NrcIsl(kIsl)
       wtot    = 0

*-->  Improve x-pos of cell by left-right info
*     and use log weighted position recon
       Do I=1,NrCells
        
        CellNr = pCellIsl(I,kIsl)
        i_Typ  = get_typ(CellNr)
*     
        Call cccXYZ(CellNr,xC,yC,zC,Lerr)
*
        xI1(kIsl) = xI1(kIsl) + ECellIsl(I,kIsl)*xC
        yI1(kIsl) = yI1(kIsl) + ECellIsl(I,kIsl)*yC
        zI1(kIsl) = zI1(kIsl) + ECellIsl(I,kIsl)*zC
*
        If (i_Typ.LE.5) emcEIsl(kIsl) = emcEIsl(kIsl) + ECellIsl(I,kIsl)

* N.T. 3/99 'max(0,...' --> max(0.,...' 
        If (i_Typ.LE.4) then
         wi  = max(0.,4.+log(ECellIsl(I,kIsl)/eIsl(kIsl)))
        Else
         wi  = max(0.,2.+log(ECellIsl(I,kIsl)/eIsl(kIsl)))
        End if
        wtot = wtot + wi         
*
        Er  = (ECellIsl(I,kIsl)+ImbCellIsl(I,kIsl))/2.0
        El  = (ECellIsl(I,kIsl)-ImbCellIsl(I,kIsl))/2.0
*
        If (Er.LT.0.0) Er = 0.0
        If (El.LT.0.0) El = 0.0
        If (El.EQ.0.0) dx = 10.
        If (Er.EQ.0.0) dx = -10.
        If (Er.NE.0.0.AND.El.NE.0.0) dx = 54./2.*log(Er/El)             
        If (Abs(dx).GT.10.0 ) dx = sign(10.0,dx)
        
************************************************************
*     Bug fix:
*
*     used POSER Nr to identify F/B/RCAL
*     this is not sensitive to shifts in the CAL position
*
*     19.5.2000 by F.Goebel
*
*       If (zC.GE.234. .OR. zC.LE.-150.) Then
*************************************************************
        if (get_fbr(CellNr).ne.2) then
*-->      FCAL or RCAL case
         xC      = xC + dx 
         xIsl(kIsl) = xIsl(kIsl) + wi*xC 
         yIsl(kIsl) = yIsl(kIsl) + wi*yC 
         zIsl(kIsl) = zIsl(kIsl) + wi*zC 
        Else
*-->      BCAL case
         d1         = Sqrt(xC**2+yC**2)
         d2(1)      = -yC/d1
         d2(2)      =  xC/d1
         xC         =  xC + dx*d2(1)
         yC         =  yC + dx*d2(2)    
         xIsl(kIsl) = xIsl(kIsl) + wi*xC
         yIsl(kIsl) = yIsl(kIsl) + wi*yC
         zIsl(kIsl) = zIsl(kIsl) + wi*zC
        EndIf
       EndDo
*     
       xI1(kIsl) = xI1(kIsl)/eIsl(kIsl)
       yI1(kIsl) = yI1(kIsl)/eIsl(kIsl)
       zI1(kIsl) = zI1(kIsl)/eIsl(kIsl)
*     
       If (wtot.GT.0) Then
        xIsl(kIsl) = xIsl(kIsl)/wtot
        yIsl(kIsl) = yIsl(kIsl)/wtot
        zIsl(kIsl) = zIsl(kIsl)/wtot
       Else
        xIsl(kIsl) = xI1(kIsl)
        yIsl(kIsl) = yI1(kIsl)
        zIsl(kIsl) = zI1(kIsl)
       EndIf
*
*-->  Find Maximum Radius of the Island on the plane perpendicular
*     to Island Ray from Vertex
       rayX = xIsl(kIsl)
       rayY = yIsl(kIsl)
       rayZ = zIsl(kIsl)-Vertex(3)

       Do I=1,NrCells

        CellNr = pCellIsl(I,kIsl)
        Call cccXYZ(CellNr,xC,yC,zC,Lerr)
        zC         = zC-Vertex(3)
*
        lenC       = Sqrt(xC**2+yC**2+zC**2)
        unit_Cx    = xC/lenC
        unit_Cy    = yC/lenC
        unit_Cz    = zC/lenC
*
        Swim       = (rayX**2+rayY**2+rayZ**2)/
     &               (unit_Cx*rayX+unit_Cy*rayY+unit_Cz*rayZ)
*
        rhoX       = Swim*unit_Cx - rayX
        rhoY       = Swim*unit_Cy - rayY
        rhoZ       = Swim*unit_Cz - rayZ
        rIsl(kIsl) = max(rIsl(kIsl),Sqrt(rhoX**2+rhoY**2+rhoZ**2))

       EndDo
     
      EndDo
*
**-->> oooooooh goody we made to the end :-)
*
      Return
      End

      Real Function probisles(dist,Itype)
*
*  Calculate the probablity that clusters should match for angle
* difference dist.  The matching for H1<->H2 is specified by itype=2
* while EMC<->H1 is specified by itype=1.
*
*  The data needed for the matching cuts
*
      Implicit none
*     
      Integer Itype
      Real ph21(4),ph1e(4),dist,slope
      Data ph21/6.527,103.7,-483.2,643.6/
      Data ph1e/11.75,100.1,-509.2,650.5/
*     
*  Get probability for EMC<->HAC1
*
      If (Itype.eq.1) then
       slope = ph1e(1)+ph1e(2)*dist+ph1e(3)*dist**2+ph1e(4)*dist**3
       probisles = exp(-slope*dist)
*
*  and now HAC1<->HAC2
*
      Elseif (Itype.eq.2) then
       slope = ph21(1)+ph21(2)*dist+ph21(3)*dist**2+ph21(4)*dist**3
       probisles = exp(-slope*dist)
      Else
       probisles = 0.
      Endif
*     
      Return
      End

