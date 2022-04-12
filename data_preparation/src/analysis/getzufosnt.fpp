c -----------------------------------------------
      subroutine getzufosnt(vtx)
c -----------------------------------------------
c ZufosNT
c 25/11/99 Write out Zufos input parameters for
c debugging purposes
C -------------------------------------------- --    
      implicit none

#include "common.inc"
#include "ZufosNT.inc"
#include "zisles.inc"

      integer l, m
      logical zfirst
      data zfirst /.true./
      real vtx(3)
      real CALPOS(3), cal_energy

      CALPOS(1)  = CAL_XP
      CALPOS(2)  = CAL_YP
      CALPOS(3)  = CAL_ZP
      cal_energy = cal_ee

c ------------------------------------------------- 
c      if (zfirst) then
c         zfirst = .false.
c         write(501,*) 'vtx1, vtx2, vtx3, NrEcells'
c         write(501,*) 'cal_xp, cal_yp, cal_zp'
c         write(501,*) 'eprob, cal_energy'
c         write(501,*) 'CellList'
c         write(501,*) '-----------------------------'
c      end if
c -------------------------------------------------
      
      do l = 1, 3
         do m = 1, 2
            zpxcal(l,m) = 0.
            zpycal(l,m) = 0.
            zpzcal(l,m) = 0.
            zufocal(l,m) = 0.
         end do
      end do
      ZEminPz = 0.
      Zgamma = 0.
      Zupt = 0.
      zuet = 0.
      ZPxBSp = 0.
      ZPyBSp = 0.
      ZPzBSp = 0.
      ZufoBSp  = 0.
      ZufoE = 0.
      ZufoPx = 0.
      ZufoPy = 0.
      ZufoPz = 0.

      call ZufosNT(vtx, NrEcells, CellList, CALPOS, EPROB,
     &     cal_energy)

c -------------------------------------------------       
c      write(501,*) cal_energy, eprob
c      write(501,*) px, py, pz, e
c      write(501,*) pxcal
c      write(501,*) pycal
c      write(501,*) pzcal
c      write(501,*) pxtrk,pytrk,pztrk,etrk
c      write(501,*) pxbsp
c      write(501,*) pybsp
c      write(501,*) pzbsp
c      write(501,*) ebsp
c      write(501,*) eminpznt,ptnt,gammant

c      write(501,*) '--------------------------------'
c      write(501,*) '--------------------------------'
c -------------------------------------------------

      do l = 1, 3
         do m = 1, 2
            ZPxCal(l, m) = PxCal(l, m)
            ZPyCal(l, m) = PyCal(l, m)
            ZPzCal(l, m) = PzCal(l, m)
            ZufoCal(l, m) = Ecal(l, m)
         end do
      end do

      ZEminPz = EminPzNT
      Zgamma = GammaNT
      ZuPt = PtNT
      zuet = Etnt
      ZPxBSp = PxBSp(2)
      ZPyBSp = PyBSp(2)
      ZPzBSp = PzBSp(2)
      ZufoBSp  = EBSp(2)
      ZufoE = E
      ZufoPx = Px
      ZufoPy = Py
      ZufoPz = Pz

      return
      end
















