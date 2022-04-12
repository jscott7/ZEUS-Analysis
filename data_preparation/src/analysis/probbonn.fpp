      Subroutine ProbBonn(vertex,NrCells,CellList,P1,Ierr)
C#*******************************************************************
C#
C@#    ProBonn : Calculate electron probability for cellclusters
C#
C#     Subroutine ProbBonn(vertex,NrCells,CellList,P1,Ierr)
C#
C#     PURPOSE: Calculate probability for Cluster Id in order
C#              to decide whether it could have been an electron
C#
C#     Communication:
C#                    vertex:   Real             vertex-information
C#                    NrCells:  Integer          Nummer of cells in
C#                                               Cluster
C#                    CellList:  Integer(NCellmax)array with
C#                                               Poser Cellnrs.
C#                    P1:       Real             Probability of being
C#                                               an elctron
C#                    Ierr:     Integer          Errorflag
C#                                               -2 NrCells <0
C#                                                2 NRCells >0
C#
C#     REMARKS: The fit has been done with the NUM4 Correction-set
C#             of the mozart T4, where the longitudinal and transvers
C#              shower-shape for electrons should be ok.
C#              A P1-value of 0.05 means, that if you cut off these
C#              condensate you will miss 5% of all real electrons
C#
C#     Author: T.Doeker
C#     Status: Preliminary
C#     Modification:
C#                      15.2.93 changed from copar
C#
C#******************************************************************
      Implicit NONE
 
#include "partap.inc"
#include "caltru.inc"
#include "elecmn.inc"
 
      Integer Nr,IC1,IC2,ICN,IFirst,NEmc,NE4Emc,Inum
      Integer Vector,I,Ierr,Nrmax,NrCells
      Real EEmc,E4Emc,x,y,z,xmax,ymax,zmax,Eges,radius,theta
      Real xabs,yabs,zabs,phi,phimax,phiabs,imbal,yabs1,Thmax
      Real P1,PE4emc,Peemc,par,EHac2,ScaEHac,ScaE4Emc
      Real vertex,Emax,Esum2,vtx
      Character*1 Cal
      Character*5 Kind
      Character*8 CName,COrder
      Logical error,cexist
      Dimension Par(4),vertex(3),vtx(3)
      Dimension Vector(5),CName(2),COrder(2)
      Data IFirst/0/
      Save IFirst,Inum
 
 
C input check
      vtx(1) = vertex(1)
      vtx(2) = vertex(2)
      vtx(3) = vertex(3)
 
      if (abs(vertex(1)).gt.10. .or.abs(vertex(2)).gt.10.
     +                       .or.abs(vertex(3)).gt.150.) then
       call VZero(vtx,3)
      endif
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C preset values
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      P1 = 0.
      Nrmax = 0
      Emax = 0.
      Ierr = 0
      EHac2 = 0.02
      Eemc = 0.
      Eges = 0.
      Esum2 = 0.
      Nemc = 0
      E4Emc = 0.
      xmax = -9999.
      ymax = -9999.
      zmax = -9999.
      imbal = 0.
      If (Ifirst.eq.0.) then
          Inum = GetInd(Caltru,'Cellnr')
          IFirst = 1
      endif
      if (NrCells.le.0) then
          Ierr = -2
          return
      endif
      if (NrCells.gt.100) then
          NrCells = 100
          Ierr = 2
      endif
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C First get Cell with highest energy
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      Do 110 I=1,NrCells
 
         if (Cexist(CellList(I))) then
 
             Caltru_Cellnr = CellList(I)
             call SelTab(Caltru,Inum,IC1,IC2)
             if(IC1.gt.IC2) then
               write(6,*) 'probonn:Cell not found. Nr.:',CellList(I)
               Ierr = 1
               goto 110
             endif
             call FetTab(Caltru,Inum,IC1)
             Call CCWhat(Caltru_CellNr,Kind,Vector,Error)
             Nr = Caltru_Cellnr
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Chimney and Hac0 -Cells are treated as Emc's
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             If ((Kind(2:2).eq.'E').or.(Kind(2:5).eq.'HAC0').or.
C!!! Beam-Pipe !!!!!!!!!!!!!!!!!!!!!!
C     .        Nr .eq. 38060  .OR.
C     .        Nr .eq. 38076  .OR.
C     .        Nr .eq. 38092  .OR.
C     .        Nr .eq. 38572  .OR.
C     .        Nr .eq. 38604  .OR.
C     .        Nr .eq. 39084  .OR.
C     .        Nr .eq. 39100  .OR.
C     .        Nr .eq. 39116  .OR.
C!!! Chimney !!!!!!!!!!!!!!!!!!!!!!!!
     .        Nr .eq. 38684  .OR.
     .        Nr .eq. 38700  .OR.
     .        Nr .eq. 38716  .OR.
     .        Nr .eq. 38732  .OR.
     .        Nr .eq. 38748  .OR.
     .        Nr .eq. 38764)       Then
              call CccXYZ(Nr,x,y,z,Error)
              call CCCaPo(x-vtx(1),y-vtx(2),z-vtx(3)
     .                                      ,radius,theta,phi)
              Phi = Phi*180/3.1416
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Position of highest energy 'Emc'-Cell
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              If (Caltru_E.gt.Emax) then
                   Emax = Caltru_E
                   Nrmax= Caltru_Cellnr
                   xmax = x
                   ymax = y
                   zmax = z
                   imbal = Caltru_imbal/Caltru_E
                   Phimax = Phi
                   Thmax  = Theta
                   E4Emc = Caltru_E
                   NE4Emc = 1
                   Phiabs = 15.
                   zabs = 7.5
                   xabs = 25.
                   yabs = 15.
                   cal = Kind(1:1)
                   if (Kind(1:1).eq.'F') yabs = 7.5
                   if ((Kind(2:5).eq.'HAC0').or.
C!!! Beam-Pipe !!!!!!!!!!!!!!!!!!!!!!
c     .                  Nr .eq. 38060  .OR.
c     .                  Nr .eq. 38076  .OR.
c     .                  Nr .eq. 38092  .OR.
c     .                  Nr .eq. 38572  .OR.
c     .                  Nr .eq. 38604  .OR.
c     .                  Nr .eq. 39084  .OR.
c     .                  Nr .eq. 39100  .OR.
C     .                  Nr .eq. 39116  .OR.
     .                  Nr .eq. 38684  .OR.
     .                  Nr .eq. 38700  .OR.
     .                  Nr .eq. 38716  .OR.
     .                  Nr .eq. 38732  .OR.
     .                  Nr .eq. 38748  .OR.
     .                  Nr .eq. 38764) yabs = 30.
              else
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C there could be an Hac0 at the side of Emax
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   yabs1 = yabs
                   if ((Kind(2:5).eq.'HAC0').or.
C!!! Beam-Pipe !!!!!!!!!!!!!!!!!!!!!!
c     .                  Nr .eq. 38060  .OR.
c     .                  Nr .eq. 38076  .OR.
c     .                  Nr .eq. 38092  .OR.
c     .                  Nr .eq. 38572  .OR.
c     .                  Nr .eq. 38604  .OR.
c     .                  Nr .eq. 39084  .OR.
c     .                  Nr .eq. 39100  .OR.
c     .                  Nr .eq. 39116  .OR.
     .                  Nr .eq. 38684  .OR.
     .                  Nr .eq. 38700  .OR.
     .                  Nr .eq. 38716  .OR.
     .                  Nr .eq. 38732  .OR.
     .                  Nr .eq. 38748  .OR.
     .                  Nr .eq. 38764)       yabs1 = 30.
 
              end if
             end if
         endif
 
  110 Continue
 
      if(nrmax.eq.0) then
*         write(6,*) 'probonn:No interesting cell found'
         ierr = 1
         return
      endif
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C now look for neighbours
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      Do 100 I=1,NrCells
 
         if (Cexist(CellList(I))) then
 
             Caltru_Cellnr = CellList(I)
             call SelTab(Caltru,Inum,IC1,IC2)
             if(IC1.gt.IC2) then
               write(6,*) 'probonn:Cell not found. Nr.:',CellList(I)
               Ierr = 1
               goto 100
             endif
             call FetTab(Caltru,Inum,IC1)
             Call CCWhat(Caltru_CellNr,Kind,Vector,Error)
             Nr = Caltru_Cellnr
             Eges = Eges + Caltru_E
             if(Kind(2:5).eq.'HAC2') Esum2 = Esum2 + Caltru_E
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Chimney and Hac0 -Cells are treated as Emc's
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             If ((Kind(2:2).eq.'E').or.(Kind(2:5).eq.'HAC0').or.
C!!! Beam-Pipe !!!!!!!!!!!!!!!!!!!!!!
c     .        Nr .eq. 38060  .OR.
c     .        Nr .eq. 38076  .OR.
c     .        Nr .eq. 38092  .OR.
c     .        Nr .eq. 38572  .OR.
c     .        Nr .eq. 38604  .OR.
c     .        Nr .eq. 39084  .OR.
c     .        Nr .eq. 39100  .OR.
c     .        Nr .eq. 39116  .OR.
C!!! Chimney !!!!!!!!!!!!!!!!!!!!!!!!
     .        Nr .eq. 38684  .OR.
     .        Nr .eq. 38700  .OR.
     .        Nr .eq. 38716  .OR.
     .        Nr .eq. 38732  .OR.
     .        Nr .eq. 38748  .OR.
     .        Nr .eq. 38764)       Then
              Eemc = Eemc + Caltru_E
              Nemc = Nemc + 1
              call CccXYZ(Nr,x,y,z,Error)
              call CCCaPo(x-vtx(1),y-vtx(2),z-vtx(3)
     .                                      ,radius,theta,phi)
              Phi = Phi*180/3.1416
              if(Nr.eq.Nrmax) goto 100
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Add Energy only if neighbour of highest EMC-Cell
C!!!!!!!!!!!!!!!!! F/RCAL then Bcal !!!!!!!!!!!!!!!!!!!!!!!
              if (((abs(xmax-x).lt.xabs)      .and.
     +             (abs(ymax-y).lt.yabs1)     .and.
     +             (Kind(1:1).ne.'B')         .and.
     +             (NE4Emc.lt.4))             .or.
     +            ((abs(Phimax-Phi).lt.phiabs).and.
     +             (abs(zmax-z).lt.zabs)      .and.
     +             (Kind(1:1).eq.'B')         .and.
     +             (NE4Emc.lt.4)))            then
                        NE4emc = NE4Emc + 1
                        E4Emc = E4Emc + Caltru_E
              end if
             end if
         end if
 100  Continue
 
C debug
C      write(6,*)'NE4EMc,E4EMC,NEMC,EEMC,Nrmax',
C     +                  NE4EMc,E4EMC,NEMC,EEMC,Nrmax
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Calculate the Scaled Parameters
C                  EHac is a function of X0 in front of
C                                        HAC-Sections
C                       and is different for Crack-events
C
C                  E4Emc depends on the Calorimeter geometry
C                                   the Imbal -Value
C                                   the Energy
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ScaEHac = 0.
      ScaE4Emc = 0.
 
      if (Eges .ne. 0.) then
        if(Cal.eq.'R'.and. cos(thmax).ne.0.)
     +         ScaEHac=(1-EEmc/Eges)*(1+6.*(-1./cos(Thmax)-1.))
        if(Cal.eq.'F'.and. cos(thmax).ne.0.)
     +         ScaEHac=(1-EEmc/Eges)*(1+6.*( 1./cos(Thmax)-1.))
        if(Cal.eq.'B'.and. sin(Thmax).ne.0.)
     +         ScaEHac=(1-EEmc/Eges)*(1+6.*(1./sin(Thmax)-1.3/25.-1.))
      end if
 
C Correction of cracks
      ScaEHac = ScaEHac*exp(-5.*abs(imbal))
C Correction of energy
      if (Eges .ne. 0.) ScaEhac = ScaEHac/sqrt(sqrt(Eges/30.))
 
C Correction of cracks
      if (EEmc .ne. 0.) then
        ScaE4Emc = (1.-E4Emc/Eemc)*exp(-abs(imbal)*3.)
 
C Geometry correction between F/BEMC and REMC
        if(Cal .ne. 'R') ScaE4Emc = ScaE4Emc/3.5
      end if
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Now use the Fit to calculat the Probability of being an electron
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C    For EEmc
         par(1) = -0.4206
         par(2) = -341.7
         par(3) = -3.557
         par(4) = -45.07
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Crack-Events have to be treated in another way
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
         if (abs(imbal).gt.0.18) then
C            For Eemc
             par(1) = -0.4571
             par(2) = -382.9
             par(3) = -1.008
             par(4) = -28.94
         end if
 
C Check (over)underflow
 
         if ((abs(par(1)+ScaEHac*par(2)) .lt. 80.) .and.
     +       (abs(par(3)+ScaEHac*par(4)) .lt. 80.)) then
                        PEEmc = exp(par(1)+ScaEHac*par(2))
     +                         +exp(par(3)+ScaEHac*par(4))
           else
                        PEEmc = 0.
         end if
 
         If (PEEmc.gt.1.) PEEmc = 1.
         If (PEEmc.lt.0.) PEEmc = 0.
 
C    For E4Emc
         par(1) = 0.2611
         par(2) = -180.4
         par(3) = -7.342
         par(4) = -165.9
 
C Check (over)underflow
 
         if ((abs(par(1)+ScaE4Emc*par(2)) .lt. 80.) .and.
     +       (abs(par(3)+ScaE4Emc*par(4)) .lt. 80.)) then
                         PE4Emc = exp(par(1)+ScaE4Emc*par(2))
     +                           +exp(par(3)+ScaE4Emc*par(4))
           else
                         PE4Emc = 0.
         end if
 
         If (PE4Emc.gt.1.) PE4Emc = 1.
         If (PE4Emc.lt.0.) PE4Emc = 0.
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Now combine these Probabilities
C it can be done because EHac and E4Emc are nearly independent
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         P1 = 0.
         if (PEEmc .le.0.000000000001) PEEmc = 0.
         if (PE4Emc.lt.0.000000000001) PE4Emc = 0.
 
         if ( PEEmc * PE4Emc .gt. 0.)
     +             P1 = PEEmc * PE4EMC * (1-log(PEEmc * PE4EMC))
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Now make a cut on the Energy in Hac2
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if(Esum2.ne.0.) then
               if(Esum2/Eges .gt. EHac2) P1 = P1*0.00001
           end if
 
Cdebug
C          write(6,*)'PEEmc,PE4emc,P1',PEEmc,PE4emc,P1
 
      Return
      End
