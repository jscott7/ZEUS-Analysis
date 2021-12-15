      logical function h1box(x,y)
      implicit none
#include "local.inc"

C#
C#  "H1" shaped box cut
C#
C#  Input: Rcal x,y electron position in cm
C#
C#  Output h1box .eq. .true. means electron passes the box cut.
C#
C#  2.12.1998 Mike Wodarczyk
C#  ------------------------------------------------------------

      real x,y

      h1box = .true.

      if ( abs(x).lt.(13+boxcutcard1)
     &.and.abs(y).lt.(7+boxcutcard1)) h1box = .false.

C
C  Feb 8,1999 F2 working group decided to widen the H1 shape to
C  the left to eliminate a region where data has more energy 
C  loss for KP than MC.  Mike Wodarczyk
C
      if ( abs(y).lt.(11+boxcutcard1)
     &.and.x.lt.(-7+boxcutcard1)
     &.and.x.gt.(-14+boxcutcard1)) h1box = .false.

C
C  new H box cut (widened on the right hand side to
C  account for data/mc beam pipe differences.)
C
      if ( abs(y).lt.(11+boxcutcard1)
     &.and.x.lt.(13+boxcutcard1)
     &.and.x.gt.(4+boxcutcard1)) h1box = .false.

      return
      end



