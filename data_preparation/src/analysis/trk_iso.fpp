C      =================================================
       SUBROUTINE TRK_ISO(trackid, cone_rad)
C      =================================================
C
c   Author: A.Quadt
c   version 1. 24.May 1998
C
C   INPUT:  trackid   -> ID of electron track in ADAMO table VCTRHL
C           delta_p   -> momentum window of tracks considered in
C                        isolation requirement
C
C   OUTPUT: cone_rad  -> cone radius of nearest non-electron track
C                        in eta-phi
C
C   ----------------------------------------------------------------

       IMPLICIT NONE

#include "partap.inc"
#include "vctrhl.inc"
#include "vcatcal.inc"
#include "const.inc"

       INTEGER trackid, I
       REAL    delta_p, cone_rad, dist
       REAL    theta_ref,  eta_ref,  phi_ref,  mom_ref
       REAL    theta_test, eta_test, phi_test, mom_test

       delta_p  = 1.0       !!! minimum p_t to be considered
       cone_rad = 9999.99
       IF (trackid.LE.0) RETURN

       CALL FETTAB(VCTRHL,  ID, trackid)
       CALL FETTAB(VCATCAL, ID, trackid)
       theta_ref = pi/2. - atan(VCTRHL_tdip)
       eta_ref   = -log(Tan(Theta_ref/2.))
       phi_ref   = atan2(VCATCAL_Py,VCATCAL_Px)
       mom_ref   = sqrt(VCATCAL_Px**2 +
     &                  VCATCAL_Py**2 +
     &                  VCATCAL_Pz**2)

       DO I=1,COUTAB(VCTRHL)
         IF (I.NE.trackid) THEN
           CALL FETTAB(VCTRHL,  ID, I)
           CALL FETTAB(VCATCAL, ID, I)
           theta_test = pi/2. - atan(VCTRHL_tdip)
           eta_test   = -log(Tan(Theta_ref/2.))
           phi_test   = atan2(VCATCAL_Py,VCATCAL_Px)
           mom_test   = sqrt(VCATCAL_Px**2 +
     &                       VCATCAL_Py**2 +
     &                       VCATCAL_Pz**2)
           IF (abs(mom_test).GT.delta_p) THEN
             dist = sqrt((eta_ref-eta_test)**2 +
     &                   (phi_ref-phi_test)**2)
             IF (dist.LT.cone_rad) cone_rad = dist
           ENDIF
         ENDIF
       ENDDO

       RETURN
       END
