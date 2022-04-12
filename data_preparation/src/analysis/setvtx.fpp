C     =======================
      SUBROUTINE SETVTX (VTX)
C     =======================
C
C     A.Quadt, August 1995
C     -----------------------------

      IMPLICIT NONE

c#include "trktup.inc"
c#include "gentup.inc"
#include "common.inc"
#include "vtxalt.inc"

      REAL    VTX(3)
      REAL    VC_CHISQ, VTX_FIX, x_fix, y_fix, z_fix
      REAL    Xerr, Yerr
      INTEGER noevt
      INTEGER TRKMIN
c
      INTEGER GOODCTD, FLAG
      REAL    BESTVTX


      WhatVertex =  0
      Z_FIX      =  0.0
      TRKMIN     =  1
      VC_CHISQ   =  5.0

      CALL GetxyVertex( RUN_NUM, X_fix, Xerr, Y_fix, Yerr, noevt )
      VTX(1) = X_FIX
      VTX(2) = Y_FIX

c      IF (VTXEL(3,3) .LT. -900.) THEN
C .......this is only satisfied before the setvtx call from sinis_st95 !
C .......set default to 0.0 cm. 
C .......(Non vertex events are biased towards FCAL)

         IF ( VCT_ZVC .GT. -900.    .AND.
     &        CHVC    .LE. VC_CHISQ .AND.
     &        NVTRKC  .GE. TRKMIN         ) THEN
            VTX(3) = VCT_ZVC
         ELSE
            VTX(3) = Z_FIX
         ENDIF

c      ELSE

C ......Multivertexing is done and VTXFLG is set: +1 if vtx found at all
C                                                    and chi2/ndf < 5 
c         GOODCTD = MOD(VTXFLG(1),2)
c         CALL BESTZVTX( VTXEL(3,3), VTXTIM, Z_FIX, GOODCTD, VTXTER,
c     +                  BESTVTX,    FLAG)
c         VTX(3)     = BESTVTX
c         WhatVertex = FLAG

c      ENDIF

C 17/09/99 TRY VTX=0 FOR THE YBUMP !!!!
C       VTX(3)     = 0

      RETURN
      END

C     ============================================================
      SUBROUTINE BESTZVTX(VTXctd, VTXcal, VTXdef, GOODctd,VTXcalE,
     +                    BESTVTX, FLAG                           )
C     =============================================================
      IMPLICIT NONE
C
      REAL     VTXctd, VTXcal, VTXdef
      INTEGER  GOODctd
      REAL     VTXcalE
      REAL     BESTVTX
      INTEGER  FLAG

      IF ( GOODCTD.EQ.1 ) THEN
         BESTVTX = VTXctd
         FLAG    = 1
      ELSEIF ( VTXcalE.LT.14. ) THEN
         BESTVTX = VTXcal
         FLAG    = 2
      ELSE
         BESTVTX = VTXdef
         FLAG    = 3
      ENDIF
      
      IF ( BESTVTX.LT.-148. .OR. BESTVTX.GT.235. ) THEN
         WRITE(*,*) 
     +        'WARNING:(setvtx) zVtx out of range (-148<z<235): ',
     +         BESTVTX,' --> put to default:',VTXdef
         BESTVTX = VTXdef
         FLAG    = 4
      ENDIF

      RETURN
      END







