      subroutine LUMIGEOM96(yes)
C     ===========================

      Implicit None

#include "common96.inc"
#include "local.inc"
C      Include 'f2isr.inc'

      Real a,b
      Logical lumigeomfirst
	integer yes

      DATA lumigeomfirst /.true./

      yes = 0

      if (lumigeomfirst.eqv..true.) then
       write (*,*) '  '
       write (*,*) '  LUMI PHOTON GEOMETRIC ACCEPTANCE CUT'
       write (*,*) '  '
       write (999,*) '  '
       write (999,*) '  LUMI PHOTON GEOMETRIC ACCEPTANCE CUT'
       write (999,*) '  '
 
       lumigeomfirst=.false.
      endif

      if (hcl_pzeg.ne.0.and.q2_tru.gt.0) then
      a = 10700.*hcl_pxeg/hcl_pzeg
Clumi_xg
      b = 10700.*hcl_pyeg/hcl_pzeg
Clumi_yg
      endif

       if (a.gt.-6.500.and.a.lt. 3.050.and.b.gt.-4.000
     & .and.b.lt. 4.000
     & .and.b+(-1.000*a).lt. 8.000
     & .and.( 0.784*a)+b.lt. 4.392
     & .and.b+(-0.976*a).gt.-4.976
     & .and.( 1.500*a)+b.gt.-12.2500) then 
         yes = 1
      endif

C
C  95 CUT     
C 
C      if(a.gt.-6.4.and.a.lt.2.85.and.b.gt.-3.7.and.b.lt.4.0.and.
C     +     (b-a).lt.7.8.and.(b-a).gt.-4.5.and.
C     +     (a+b).lt.4.5.and.(a+b).gt.-9.0) then
C         lumigeom = .true.
C      endif
C
C       if(a.gt.-6.75.and.a.lt.3.85.and.b.gt.-4.0.and.b.lt.4.0.and.
C     +     (b-1.145*a).lt.8.58.and.(b-0.933*a).gt.-4.79.and.
C     +     (0.792*a+b).lt.4.0.and.(1.05*a+b).gt.-8.988) then
C         lumigeom = .true.
C       endif
C
C       if(a.gt.-5.75.and.a.lt.2.85.and.b.gt.-3.0.and.b.lt.3.0.and.
C     +     (b-1.229*a).lt.7.916.and.(b-0.933*a).gt.-3.79.and.
C     +     (0.7192*a+b).lt.3.0.and.(1.05*a+b).gt.-7.988) then
C         lumigeom = .true.
C       endif
C
C       if (a.eq.0.0.or.b.eq.0.0) then
C         lumigeom = .false.
C       endif

      
      return
      end
