C     ===========================================
      SUBROUTINE MY_GPXYZ(CPX,CPY,CPZ,CET,RAD,XV)
C     ===========================================
C     Auxiliary routine for angle calculation.
C
C     A.Quadt, May 1997.
C     -----------------------------------------

      IMPLICIT NONE

#include "caltru.inc"
#include "partap.inc"

C
C -------> Caltru is Fettab'd
C
      INTEGER    IERR
      REAL       ZP,YP,CPX,CPY,CPZ,CET,RAD,XV(3),THETA,PHI,
     1           SINPHI,SINTH,COSTH,COSPHI,XP

      CALL CCCXYZ (CALTRU_CELLNR,XP,YP,ZP,IERR)
      CALL CCCAPO (XP-XV(1),YP-XV(2),ZP-XV(3),RAD,THETA,PHI)
      RAD    = SQRT((XP-XV(1))**2+(YP-XV(2))**2+(ZP-XV(3))**2)
      SINPHI = SIN(PHI)
      COSPHI = COS(PHI)
      SINTH  = SIN(THETA)
      COSTH  = COS(THETA)

      CPX    = CALTRU_E*SINTH*COSPHI
      CPY    = CALTRU_E*SINTH*SINPHI
      CPZ    = CALTRU_E*COSTH
      CET    = CALTRU_E*SINTH

      RETURN
      END

