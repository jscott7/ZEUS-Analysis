*      ==================
       SUBROUTINE MULOOP(rver,  igmu)
*      ==================
       IMPLICIT NONE
* -----------------------------
* Inputs-
*         rver(3)             -  The primary vertex x,y,z position
*         common block MUOTUP -  Filled by previous CALL to MUOFIL
* Output-
*         igmu = 1 -> muon candidate found
*              = 0 -> no muon candidate found
* -----------------------------
#include "muotup.inc"
       INTEGER igmu
       REAL    RVER(3)
*
       REAL VTXCA,BMUoppCA(2)
*
       REAL BMUCAVTXMAX,BMUOPPCAMIN,BMUSEPMAX
       PARAMETER  ( BMUCAVTXMAX = 100. )
       PARAMETER  ( BMUOPPCAMIN = 50.  )
       PARAMETER  ( BMUSEPMAX   = 200. )
c
       REAL    Rseg1(3),RDCOS1(3),RCA(3),POINT_SEP
       REAL    RDCOS2(3),RSEG2(3),RCA1(3),RCA2(3)
       REAL    DCA1to2,DCA2to1
       REAL    vtxDCANEW,combinedsepmin,combinedsep
       INTEGER i,j
       LOGICAL FOUND_OPPOSITES
c
       LOGICAL first
       DATA first /.true./
       SAVE first
c ---------
       IF (first) THEN
         first=.false.
         write(6,*)' --------- GAZZAMU -----------'
         write(6,*)' PARAMETER LIST: 24/6/97'
         write(6,*)' BMUCAVTXMAX ',BMUCAVTXMAX
         write(6,*)' BMUOPPCAMIN ',BMUOPPCAMIN
         write(6,*)' BMUSEPMAX   ',BMUSEPMAX
         write(6,*)' ---------------------------'
       ENDIF
c ---------
c Initiate initial initializations.
c
       igmu            =     0
       igmuon          =     0
       vtxCA           =    -1.
       BMUoppCA(1)     = 10000.
       BMUoppCA(2)     = 10000.
       combinedsepmin  = 10000.
       FOUND_OPPOSITES = .false.
c
c Loop over all bmuon tracks.
c
       DO i=1,nbmtk
c
c Copy X,Y,Z info to R(3) vectors to ease  manipulation.
c
           Rseg1(1)  = Xseg(i)
           Rseg1(2)  = Yseg(i)
           Rseg1(3)  = Zseg(i)
c
           RDCOS1(1) = XDcos(i)
           RDCOS1(2) = YDcos(i)
           RDCOS1(3) = ZDcos(i)
c
c Find perigee to vertex for those with inner & outer bmu hits only
c Store the largest DCA to vertex of all the tracks in VTXCA.
c
         IF (INVp(i).lt.1000.) THEN
           CALL BMUPERI(RVER,Rseg1,RDcos1,RCA)
           vtxDCANEW = POINT_SEP(rver,RCA)
*
           if ( vtxDCANEW.gt.vtxCA ) vtxCA = vtxDCANEW
         ENDIF
c
c For each track loop over all others.  Check for separated
c tracks pointing to each other.  Classic cosmic, man.
c
        DO j=1,nbmtk
c
c Copy X,Y,Z info to R(3) vectors to ease  manipulation.
c
          Rseg2(1)  = Xseg(j)
          Rseg2(2)  = Yseg(j)
          Rseg2(3)  = Zseg(j)
*
          RDCOS2(1) = XDcos(j)
          RDCOS2(2) = YDcos(j)
          RDCOS2(3) = ZDcos(j)
*
*  skip same-track comparison and tracks too close together.
*
          if (j.eq.i .or.
     +     POINT_SEP(Rseg1,Rseg2).lt.BMUSEPMAX) goto 10
*
           CALL BMUPERI(Rseg2,Rseg1,RDcos1,RCA1)
           DCA1to2 = POINT_SEP(Rseg2,RCA1)
*
           CALL BMUPERI(Rseg1,Rseg2,RDcos2,RCA2)
           DCA2to1 = POINT_SEP(Rseg1,RCA2)
*
           combinedsep = DCA1to2 + DCA2to1
           IF   ( DCA1to2.lt.bmuoppcamin
     +       .and.DCA2to1.lt.bmuoppcamin
     +       .and.combinedsep.lt.combinedsepmin) THEN
*
             FOUND_OPPOSITES = .true.
*
             combinedsepmin  = combinedsep
             BMUoppCA(1)     = DCA1to2
             BMUoppCA(2)     = DCA2to1
           ENDIF
c --
 10     CONTINUE
       ENDDO
c
       ENDDO
c
c       BMUOK=(.not.FOUND_OPPOSITES).and.VTXCA.lt.BMUCAVTXMAX
c
       IF (FOUND_OPPOSITES .OR. VTXCA.gt.BMUCAVTXMAX) THEN
        igmu=1
       ELSE
        igmu=0
       ENDIF
c
c Set duplicate flag in common block
c
       IGMUON = IGMU
c
       END
 
 
c ----------
c      =======================
       REAL FUNCTION POINT_SEP(X,Y)
c      ======================
       IMPLICIT NONE
c
c Auxilliary function of GAZZAMU.  Called by MULOOP routine.
c Returns distance between two points. Not clever, but effective.
c
       REAL X(3),Y(3)
       POINT_SEP=sqrt((X(1)-Y(1))**2+(X(2)-Y(2))**2+(X(3)-Y(3))**2)
       END
c ----------
c      ==================
       SUBROUTINE BMUPERI(coord,Rseg,Rdcos,RCA)
c      ==================
       IMPLICIT NONE
c
c Auxilliary subroutine of GAZZAMU. Called by MULOOP routine.
c
c Finds closest approach of BMUON tracks to the
c given coords (e.g. the vertex
c                    or a high Et island
c                    or another muon track segment).
c
c Inputs.  Coord(3) R : coordinate to find closest approach to.
c          Rseg(3)  R : X,Y,Z of BMU track segments
c          Rdcos(3) R ; direction cosines of Bmu track.
c
c Output.  RCA(3)   R : coordinate of point of closest appraoch of
c                           extrapolated track
c
       REAL Coord(3),tzero,Rseg(3),RDcos(3),RCA(3)
c
           tzero =  ( (coord(1)-Rseg(1))*RDCOS(1)
     &              + (coord(2)-Rseg(2))*RDCOS(2)
     &              + (coord(3)-Rseg(3))*RDCOS(3) )
c
          RCA(1) = Rseg(1) + tzero*RDcos(1)
          RCA(2) = Rseg(2) + tzero*RDcos(2)
          RCA(3) = Rseg(3) + tzero*RDcos(3)
       END
 
