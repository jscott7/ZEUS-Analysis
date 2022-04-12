C     =================
      SUBROUTINE TRKFIL
C     =================
C     This routine reconstructs the event tracking vertex,
C     the number of tracks, the number of vertex tracks and the
C     chi**2/ndf for the vertex fit. At present the vertex
C     parameters are stored for `CTD-only' and `regular mode',
C     all other track parameters (cal track matching etc.) are
C     determined from `regular tracking' mode.
C
C     A.Quadt, May  1997.
C              July 1997, CTD and global tracking modes added.
C     O.Ruske, January 1998.
C              added the three angles needed to determine the
C              vertex distribution: theta_perpf,b,r
C              (from CTD only)  see note book II, page 81
C     ----------------------------------------------------------

      IMPLICIT NONE

#include "partap.inc"
#include "vctvtx.inc"
#include "vctpar.inc"
#include "vctrhl.inc"
#include "zdskey.inc"
#include "zrevt.inc"
#include "vcevtctd.inc"
#include "vcevtreg.inc"
C
c#include "zvttup.inc"
#include "common.inc"
C
      INTEGER     NO_USE, Ierr, NrGoodCells, VzErr, I
      DATA        NO_USE / -999. /
      REAL        COMPF, COMPB, COMPR, PI
      REAL        VtxtFCAL
      REAL        VtxtFTD(3)
 
      PI =  ACOS(-1.0)
C ---------------------------------------------
C --- Unpack the tracking tables (CTD only) ---
C ---------------------------------------------

      IF (CouTab(VCEVTCTD).ne.0) CALL VCGETCTD(Ierr)

C -----------------------------------------
C --- reconstructed vertex from VCTRACK ---
C -----------------------------------------
      CHVCC   = NO_USE
      VCT_XVC = NO_USE
      VCT_YVC = NO_USE
      VCT_ZVC = NO_USE 
      NVTRKC  = COUTAB(VCTPAR)
      NTRKC   = COUTAB(VCTRHL)
       
      IF (COUTAB(VCTVTX).GT.0) THEN
            CALL FETTAB(VCTVTX,ID,coutab(vctvtx))
            VCT_XVC = VCTVTX_V(1)
            VCT_YVC = VCTVTX_V(2)
            VCT_ZVC = VCTVTX_V(3)
            IF (VCTVTX_NDF .GT.0.AND.
     &          VCTVTX_CHI2.GT.0) THEN
              CHVCC = VCTVTX_CHI2/FLOAT(VCTVTX_NDF)
            ENDIF
        ENDIF    
c      ENDIF

C ------------------------------------------------------------
C --- Unpack the tracking tables (regular mode; incl. FTD) ---
C ------------------------------------------------------------
c      IF (CouTab(VCEVTREG).ne.0) CALL VCGETREG(Ierr)

C -----------------------------------------
C --- reconstructed vertex from VCTRACK ---
C -----------------------------------------
c      CHVC   = NO_USE
c      VCT_XV = NO_USE
c      VCT_YV = NO_USE
c      VCT_ZV = NO_USE
c      IF (IErr.EQ.0) THEN
c        NVTRK  = COUTAB(VCTPAR)
c        NTRK   = COUTAB(VCTRHL)

c        IF (COUTAB(VCTVTX).GT.0) THEN
c            CALL FETTAB(VCTVTX,ID,1)
c            VCT_XV = VCTVTX_V(1)
c            VCT_YV = VCTVTX_V(2)
c            VCT_ZV = VCTVTX_V(3)
c            IF (VCTVTX_NDF .GT.0.AND.
c     &          VCTVTX_CHI2.GT.0) THEN
c              CHVC = VCTVTX_CHI2/FLOAT(VCTVTX_NDF)
c            ENDIF
c        ENDIF
c      ENDIF


C ----------------------------------------
C --- FCAL timing vertex from R.Pawlak ---
C ----------------------------------------
      FCAL_VTX  = VtxtFCAL(NrGoodCells,VzErr)
      FCAL_VTXE = VzErr

C ----------------------------------------
C --- FTD vertex from M.Eckert
C ----------------------------------------
      FTD_VTX  = -999.9
      FTD_VSEG = 0
      CALL FTDVTX(VTXTFTD,FTD_VSEG)
      FTD_VTX = VtxtFTD(3)


      RETURN
      END


