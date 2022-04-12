
C     =================
      SUBROUTINE BEAMFIL
C     =================
C     This routine fills general event quantities such as
C     run and event number, the year the data was taken
c     (or simulated) and the beam energies.
C
C     A.Quadt, May 1997.
C     ----------------------------------------------------

      IMPLICIT NONE

#include "partap.inc"
#include "zrevt.inc"
#include "zdskey.inc"
#include "fmckin.inc"
#include "o4sbor.inc"
c#include "ctrlcm.inc"
#include "common.inc"

      INTEGER EVT_COUNTER, I
      SAVE    EVT_COUNTER
      DATA    EVT_COUNTER / 0 /
      REAL    DEF_EE,DEF_EP

C ---------------------
C --- Beam energies ---
C ---------------------
      IF (YEAR.EQ.1993) THEN
          DEF_EE = 26.667
          DEF_EP = 819.92
      ELSEIF (YEAR.EQ.1994) THEN
          DEF_EE = 27.584
          DEF_EP = 819.901
      ELSEIF (YEAR.EQ.1995) THEN
          DEF_EE = 27.584
          DEF_EP = 819.901
      ELSEIF (YEAR.EQ.1996) THEN
          DEF_EE = 27.584
          DEF_EP = 819.901
      ELSEIF (YEAR.EQ.1997) THEN
          DEF_EE = 27.584
          DEF_EP = 819.901
      ENDIF

      IF (.not.montecarlo) THEN

        IF (COUTAB(O4SBOR).GT.0) THEN
          CALL FETTAB(O4SBOR,ID,1)
          EE     = O4SBOR_Momnz
          EP     = O4SBOR_Mompz
          DEF_EE = EE
          DEF_EP = EP
        ELSE
          EE     = DEF_EE
          EP     = DEF_EP
        ENDIF

      ELSE

        IF (COUTAB(FMCKIN).GT.0) THEN
          CALL FETTAB(FMCKIN,ID,1)     ! fetch initial electron energy
          EE     = FMCKIN_P(4)
          CALL FETTAB(FMCKIN,ID,2)     ! fetch initial proton energy
          EP     = FMCKIN_P(4)
          DEF_EE = EE
          DEF_EP = EP
        ELSE
          EE     = DEF_EE
          EP     = DEF_EP
        ENDIF

      ENDIF



      RETURN
      END


