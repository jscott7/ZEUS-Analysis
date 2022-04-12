C     =================
      Subroutine VtxAlt
C     =================
      Implicit None
C
C Get alternative vertices:
C     FCAL Timing Vertex and FDET vertex.
C
C
#include "vtxalt.inc"
c
      Real    VtxtFCAL
      Real    VzErr
      Integer NrGoodCells

      Real    Vtx(3)
      Integer mlt

C--
C Radek's FCAL TIMING VERTEX:
C --
      VtxTim = VtxtFCAL(NrGoodCells,VzErr)
      VtxTEr = VzErr

C--
C FTD VERTEX:
C--
      Call FTDVTX( Vtx, mlt)
      VtxFtd = Vtx(3)


      Return
      End
