C     ======================================
      subroutine QEDComp( vtx, IFlag, esum )
C     ======================================
C
C------------------------------------------------------------------
C
C Purpose: Find QED Compton events
C
C Input: vtx real(3): event vertex
C
C Output: IFlag integer: 0 = no QED compton
C                        1 = QED compton
C         Esum  real:    Energy of electron and gamma
C
C Author: T.Doeker
C
C Date: 17.08.93
C
C Status: Test
C
C----------------------------------------------------------------
 
      implicit none
 
#include "cconsa.inc"
#include "partap.inc"
#include "zdskey.inc"
#include "tctrak.inc"
#include "vctrhl.inc"
#include "elecmn.inc"
 
      real    vtx(3), esum
      real    Emin, Esummin, Esummax, Esmall,prob
 
      integer N,I,IFlag, NTc, NVc,ierr,NCell
 
      parameter ( Emin = 1., Esummin = 15.5, Esummax = 30. )
 
C...preset values
 
      NTc = 0
      NVc = 0
      N = 0
      IFlag = 0
      Esum  = 0.
      Esmall = 0.
 
C...loop over condensates
 
      do I=1,Coutab(CConSa)
         call FetTab(CConSa,Id,I)
         call concells(CConSa_Id,NCell,CellList,Ierr)
         call probbonn(vtx,NCell,CellList,Prob,Ierr)
         if (CConSa_E.gt.Emin            .and.
c     +        CConSa_Class.eq.'electron' .and.
     +        prob .ge. 0.01              .and.
     +        (CConSa_z.lt.230. .or.
     +         abs(CConSa_x).gt.30. .or. abs(CConSa_y).gt.30.)) then
             Esum = Esum + CConSa_E
             N = N + 1
             if (N.eq.1) then
             endif
             if (N.eq.2) then
             endif
 
         else
             Esmall = Esmall + CConSa_E
         endif
      enddo
 
 
C...loop over tracks
 
      do I=1,Coutab(Tctrak)
         call FetTab(TcTrak,Id,I)
         if (TcTrak_NDofF .ge. 7) then
            if (TcTrak_chi2/TcTrak_NDofF .ge.4.) Ntc = Ntc + 1
         endif
      enddo
      do I=1,Coutab(vctrhl)
         call FetTab(VcTrhl,Id,I)
         if (Vctrhl_Ndf .ge. 7) then
            if (VcTrhl_chi2/VCtRhl_Ndf .lt.4.) Nvc = Nvc + 1
         endif
      enddo
 
C...select
 
      if (Esmall.gt.Emin)                                return
      if (N .ne. 2)                                      return
      if (Esum.gt.Esummax .or. Esum.lt.Esummin)          return
      if (Nvc.gt.1)                                      return
 
      IFlag = 1
C...debug
      print *, '    * QED compton *', Esmall,Esum,N,coutab(tctrak)
C...debug
      return
      end
