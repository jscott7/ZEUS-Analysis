      subroutine z_cIsles(vtx,z_Ncells,z_cE,z_cPNr,z_cImb,z_nI,CorT,
     &                     ierr)
       IMPLICIT NONE
*-------------------------------------------------------------------------
*
*      INPUT:
*               vtx(3)       <--> Reconstructed Vertex
*               z_Ncells     <--> number  of cells 
*               z_cE(*)      <--> Energy  of cell  
*               z_cPNr(*)    <--> PoserNr of cell  
*               z_cImb(*)    <--> Imbl    of cell  
*               CorT         <--> cell(1,2,11,12) or tower(20) Islands
*
*     OUTPUT:   z_nI         <--> number of islands   
*       
*       Ierr:   1  <--> No islands found
*               0  <--> Everything is cosher 
*              -1  <--> nIslands > nI_Max
*              -10 <--> Problem in cellIsland
*              
*  Special thanks to Bruce Straub for many useful and 
*  original discussions about cellIslands and clustering in general.
*
*  Tower Islands were developed by Larry Wai and could be obtained 
*  by calling this routine with CorT flag = 20.
*
*  Cell Islands:
*  CorT mode =1 cells that match at the corner are YES connected
*  CorT mode =2 cells that match at the corner are NOT connected
*
*  CorT mode =11 cells that match at the corner are YES connected,
*                no connection in depth (2-D mode)
*  CorT mode =12 cells that match at the corner are NOT connected, 
*                no cennection in depth (2-D mode)
*
*  Author:   G. Briskin 24.07.96 (Tel Aviv University)    
*
*__________________________________________________________________________
*
#include "zisles.inc"
!  /*   ! Island  */
*
*==>  Input variables
      Real    vtx(3)
      Integer z_Ncells , z_cPNr(*)
      Real    z_cE(*)  , z_cImb(*)
      Integer CorT                 !1,2,11,12=cellIslands, 20=towerIslands

*==>  Output variables
      Integer z_nI,Ierr

*==>  Internal Island variables
      Integer cPNr( Nr_cMax),cIsl( Nr_cMax)
      REAL    cE(   Nr_cMax),cImb( Nr_cMax)

*==>  cellIsland clustering mode
      Integer CorTIsl,Imode

*==>  sorting index
      Integer INDEX(Nr_cMax),RevInd(Nr_cMax)

*==>  Counters
      INTEGER kIsl,kCell,I,NrCells
      
*==>  CCrecon things
      Integer CellNr,icel
      Logical Lerr

*==>  Position Improvement
      Real    xC, yC, zC
      Real    dx, Er, El, d1, d2(2), wi, wtot

*==>  Variables for projection on perpendicular plane
      Real    rayX,rayY,rayZ
      Real    lenC,unit_Cx,unit_Cy,unit_Cz
      Real    rhoX,rhoY,rhoZ
      Real    Swim

*==>  function declarations
      Integer cellIsland  
      Integer get_fbr, get_mod, get_tow, get_typ
*
      get_fbr(CellNr) = ISHFT(CellNr,-14)+1
      get_Mod(CellNr) = ISHFT(IAND(CellNr,x'3E00'),-9)+1
      get_Tow(CellNr) = ISHFT(IAND(CellNr,x'1F0' ),-4)+1
      get_typ(CellNr) = ISHFT(IAND(CellNr,x'E'   ),-1)
*
*==> initialize
      Ierr = 0
      z_nI = 0
      nIsl = 0

*==>  First check if we are running in an allowed mode
      If    (CorT.EQ.1) Then
       CorTIsl  = 1
       Imode    = 1
      ElseIf(CorT.EQ.2) Then
       CorTIsl  = 1
       Imode    = 2
      ElseIf(CorT.EQ.11) Then
       CorTIsl  = 1
       Imode    = 11
      ElseIf(CorT.EQ.12) Then
       CorTIsl  = 1
       Imode    = 12
      ElseIf(CorT.EQ.20) Then
       CorTIsl  = 2       
      Else
       STOP 'z_cIsles: Wrong Processig Mode !!!'
      EndIf

*==>  Zero User Info
      Do I=1,nI_Max
       eIsl(    I) = 0.
       xIsl(    I) = 0.
       yIsl(    I) = 0.
       zIsl(    I) = 0.
       rIsl(    I) = 0.
       xI1(     I) = 0.
       yI1(     I) = 0.
       zI1(     I) = 0.
       emcEIsl( I) = 0.
       NrcIsl(  I) = 0
       IslTyp(  I) = 0
      EndDo

*==>  Prepare input to cellisland(i.e. desending in energy)     
      If     (CorTIsl.EQ.1) Then

       CALL SORTZV(z_cE,INDEX,z_Ncells,1,1,0)
       Do I=1,z_Ncells
        cPNr(  I) = z_cPNr(INDEX(I))
        cE(    I) = z_cE(  INDEX(I))
        cImb(  I) = z_cImb(INDEX(I))
        RevInd(I) =        INDEX(I)
       EndDo

      ElseIf (CorTIsl.EQ.2) Then

       Do I=1,z_Ncells
        cPNr(I)      = z_cPNr(I)
        cE(  I)      = z_cE(  I)
        cImb(I)      = z_cImb(I)
        RevInd(I) =           I
       EndDo

      EndIf

*==>  Construct Islands..........................................
      If     (CorTIsl.EQ.1) Then

       z_nI = cellIsland(Imode, z_Ncells, cPNr, cIsl)
       If (z_nI.LT.0) Then
        z_nI =  0
        Ierr = -10
        Return
       End If
       
      ElseIf (CorTIsl.EQ.2) Then
       
       Call Islands(z_nI,cIsl,cPNr,cE,z_Ncells,Ierr)
       If ( Ierr.NE.0 ) then
        z_nI =  0
        Ierr = -10
        Return
       End If
       
      EndIf
*
      If (z_nI.EQ.0) Then
       Ierr = 1
       Return
      ElseIf (z_nI.GT.nI_Max) Then
       z_nI =  0
       Ierr = -1
       Return
      EndIf

*==>  Store Number of cellIslands 
      nIsl = z_nI

*==>> find island energy and position
      Do kCell=1,z_Ncells 
       kIsl                          = cIsl( kCell)

       NrcIsl( kIsl)                 = NrcIsl(kIsl) + 1
       zIsl_Cell(NrcIsl(kIsl),kIsl)  = RevInd(kCell)
       
       eIsl(                   kIsl) = eIsl(kIsl) + cE(kCell)            
       pCellIsl(  NrcIsl(kIsl),kIsl) = cPNr(kCell)
       ECellIsl(  NrcIsl(kIsl),kIsl) = cE(  kCell) 
       ImbCellIsl(NrcIsl(kIsl),kIsl) = cImb(kCell) 
      EndDo

*==>  get Island position...............................................
      Do kIsl=1,z_nI
*
       NrCells      = NrcIsl(kIsl)
       CellNr       = pCellIsl(1,kIsl)
       icel         = get_typ(CellNr)

*==>   This is needed only for coneIslands (it is not used in any other place)
       If     (icel.LE.5) Then
        IslTyp(kIsl) = 1000
       ElseIf (icel.EQ.6) Then 
        IslTyp(kIsl) = 2000
       ElseIf (icel.EQ.7) Then 
        IslTyp(kIsl) = 3000
       EndIf

       wtot         = 0

*==>   Improve x-position of cell by left-right info + log weighting......
       Do I=1,NrCells

        CellNr = pCellIsl(I,kIsl)
        icel   = get_typ(CellNr)
*
        CALL cccXYZ(CellNr,xC,yC,zC,Lerr)
*         
        xI1(kIsl) = xI1(kIsl) + ECellIsl(I,kIsl)*xC
        yI1(kIsl) = yI1(kIsl) + ECellIsl(I,kIsl)*yC
        zI1(kIsl) = zI1(kIsl) + ECellIsl(I,kIsl)*zC
*
        If (icel.LE.5) emcEIsl(kIsl)=emcEIsl(kIsl)+ECellIsl(I,kIsl)
* 
* N.T. 3/99 'max(0,...' --> max(0.,...' 
        If (icel.LE.4) then
         wi  = max(0.,4.+log(ECellIsl(I,kIsl)/eIsl(kIsl)))
        Else
         wi  = max(0.,2.+log(ECellIsl(I,kIsl)/eIsl(kIsl)))
        End if
        wtot = wtot + wi
*        
        Er   = (ECellIsl(I,kIsl)+ImbCellIsl(I,kIsl))/2.0
        El   = (ECellIsl(I,kIsl)-ImbCellIsl(I,kIsl))/2.0
*
        If (Er.LT.0.0) Er = 0.0
        If (El.LT.0.0) El = 0.0
        If (El.EQ.0.0) dx = 10.
        If (Er.EQ.0.0) dx = -10.
*!!!! The lambda=54. is from Kruger DESY report, 
*!!!! The lambdas are different and their values should be looked up in the ZEUS-Note 97-019
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
*==>     FCAL or RCAL case
         xC      = xC + dx 
         xIsl(kIsl) = xIsl(kIsl) + wi*xC 
         yIsl(kIsl) = yIsl(kIsl) + wi*yC 
         zIsl(kIsl) = zIsl(kIsl) + wi*zC 
        Else
*==>     BCAL case
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
       xI1( kIsl) = xI1(kIsl)/eIsl(kIsl)
       yI1( kIsl) = yI1(kIsl)/eIsl(kIsl)
       zI1( kIsl) = zI1(kIsl)/eIsl(kIsl)
*
       If (wtot.GT.0.AND.NrCells.LE.30) Then
        xIsl(kIsl) = xIsl(kIsl)/wtot
        yIsl(kIsl) = yIsl(kIsl)/wtot
        zIsl(kIsl) = zIsl(kIsl)/wtot
       Else
        xIsl(kIsl) = xI1(kIsl)
        yIsl(kIsl) = yI1(kIsl)
        zIsl(kIsl) = zI1(kIsl)
       EndIf
*  
*==>  Find Maximum Radius of the Island on the plane perpendicular
*     to Island Ray from Vertex
       rayX = xIsl(kIsl)
       rayY = yIsl(kIsl)
       rayZ = zIsl(kIsl)-vtx(3)

       Do I=1,NrCells

        CellNr = pCellIsl(I,kIsl)
        Call cccXYZ(CellNr,xC,yC,zC,Lerr)
        zC         = zC-vtx(3)
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
*===> oooooooh goody we made to the end :-)
*
      End

