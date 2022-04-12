      subroutine ConCells(Cid,NrCells,CellList,Ierr)
C#**********************************************************
C#
C#      subroutine ConCells(Cid,NrCellsCellList,Ierr)
C#
C@#  routine returns the Cellnr. for a  Condensat
C#
C#   this routine fills the Cellnrs of a Condensat into
C#   the CellList array
C#
C#     called by user
C#
C#   Input:  Cid       Condensat Id
C#
C#   Output: CellList  Integer(NCellmax)
C#           NrCells   Nr. of Cells in Condensate
C#           Ierr      Integer      0 o.k.
C#                                 -1 no condensates found
C#                                  1 more then Ncellmax Cells
C#
C#  Author: Thomas Doeker (vxdesy::doeker)
C#  Date:   30.4.93
C#  Status: test
C#
C#*********************************************************
 
      Implicit None
 
#include "caltru.inc"
#include "cconsa.inc"
#include "partap.inc"
#include "elecmn.inc"
 
      Integer CId,Ierr,NrCells,
     +        J,Inum,IFirst
      Data Ifirst/0/
      Save Ifirst,INum
 
       If (IFirst.eq.0.) then
           Inum = GetInd(Caltru,'E')
           IFirst = 1
       endif
C     preset values
 
      Ierr   = 0
      NrCells= 0
 
      if (CId.le.0) then
          Ierr = -1
          return
      endif
 
      call Vzero(CellList,NCellMax)
 
C     select condensate
 
      do 110 J=1,CouTab(Caltru)
 
C     first take Condensates with maximum Energy
 
          call FetTab(Caltru,Inum,Coutab(Caltru)+1-J)
 
C     select the cells
 
         If (Caltru_CConSa.eq. CId) then
             If (NrCells.lt.NCellmax) then
                 NrCells = NrCells + 1
                 CellList(NrCells) = Caltru_CellNr
             else
                 Ierr = 1
             endif
         endif
 
  110 continue
 
      return
      end
