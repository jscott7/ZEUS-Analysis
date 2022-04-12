      subroutine preptracks(ierr) 
*----------------------------------------------------------------------
*     PURPOSE:
*     ========
*     calls vceaze in multi-vertex mode 
*     checks if any vceaze error has occured, if yes, errorflag 
*     gets set to 1, otherwise ierr = 0
*
*----------------------------------------------------------------------

      implicit none

#include "partap.inc"
#include "vcevtctd.inc"

      integer ierr
      integer ivcerr

      integer jinit, jtrack, jfit, jdedx, jvertex, jerr
      logical Lfirst
      data Lfirst  /.true./

      ierr = 0   ! everthing is ok
c
      if (Lfirst) then
c
*     --- VCEAZE setup
         jinit   = 0
         jtrack  = 0
         jfit    = 0
         jdedx   = 0
         jvertex = 1
         Lfirst = .false.
      endif
c
      call vceaze(jinit, jtrack, jfit, jdedx, jvertex, jerr)
c     
      if (coutab(vcevtctd).gt.0) then
         call vcgetctd(ivcerr)
         if (ivcerr.ne.0)  then 
            ierr = 1
         endif
      endif

      return
      end



