c ---------------------------------------------
      subroutine getpres(vtx)
c ---------------------------------------------
C PRESAMPLER information 
C --------------------------------
      implicit none

#include "common.inc"

      real vtx(3)
      real CALPOS(3)
      integer Ierr
      integer CENTRAL_TILE
      REAL PSAME(3)
      
      CALPOS(1) = CAL_XP
      CALPOS(2) = CAL_YP
      CALPOS(3) = CAL_ZP

      PRSE1 = 0.
      PRSE2 = 0.
      PRSE3 = 0.

      call PRCALIB(Ierr)
      CALL PRCLUS(CALPOS, vtx(3), PSAME, CENTRAL_TILE, Ierr)
      IF (Ierr .GE. 0) THEN
         PRSE1 = PSAME(1)
c --------------------------------------------------------
c PSNAME(2) -- presampler energy, use later for correction
c --------------------------------------------------------
         PRSE2  = PSAME(2) 
c--------------------------------------------------------
         PRSE3  = PSAME(3)
      ENDIF

      return
      end
