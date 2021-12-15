      subroutine best_theta(best,cal,srtd,hes,tdip,tvctpar,vtx,mode
     &     ,montecarlo)
C
C best_theta returns the best reconstructed angle for the scattered electron
C given the CAL, SRTD, HES and track dip angle for the electron and the 
C measured vertex.
C
C INPUT:
C    Real cal(4)  
C         cal(1) = x cal from elecpo
C         cal(2) = y cal from elecpo
C         cal(3) = z cal from elecpo
C         cal(4) = electron energy (not used by best_theta)
C
C    Real srtd(4)  
C         srtd(1) = x from srtdelec
C         srtd(2) = y from srtdelec
C         srtd(3) = z from srtdelec (not used)
C         srtd(4) = srtd mips 
C                   (if mips is in the range 0.4 to 1000 then srtd can be used)
C                   (NO srtd mips are returned by srtdelec if less than 0.4)
C
C    Real hes(4)  
C         hes(1) = x from hesclu
C         hes(2) = y from hesclu
C         hes(3) = z from hesclu (not used)
C         hes(4) = hes mips 
C                   (if mips is in the range 5 to 1000 then hes can be used)
C                   (NO HES mips are returned by hesclu if less than 5 mips)
C
C    Real tdip   = tracking dip value from VCTRHL
C
C    Real tvctpar  = tracking angle  from VCTPAR (not used yet set to -1. )
C
C    Real vtx(3) = vertex (x,y,z) position
C
C    Integer mode = mode to run best_theta
C             mode = 1 means just calculate the best angle and quit
C                    if called twice for the same event, the calculation
C                    is not re-done, but the old result is returned for speed.
C                    (Recommended)
C          
C             mode = 2 means Do full calculation each time.  Calculate
C                    all angles (cal,hes,tdip...) and leave them in
C                     common block.  If called more than once per event,
C                     recalculate each time.
C
C    Logical Montecarlo = True if montecarlo
C
C
C Author: Mike Wodarczyk  13-June-1998   Version 1.0
C
C Revision 1:  Mike Wodarczyk 15-June-1998 Version 1.1
C     Modified HES crack code so that Module 12 works.
C     Modified HES Z in data and MC.
C     Modified SRTD outer edge to handle SRTD crack along 
C              Rcal half edges.
C
C Revision 2: Mike Wodarczyk 18-June-1998 Version 1.2
C     Fixed Lower y srtd edge cut.
C
C Revision 3: Jonathan Scott 21-April-1999
C     Require tdip input to be tracking angle (not cot(theta))
C
      Implicit None
#include "best_theta.inc"
      Real cal(4),srtd(4),hes(4),tdip,tvctpar,vtx(3)
      Integer mode
      Logical Montecarlo
      Real best

      real dataRCALz,mcRCALz, rcalz
      parameter ( dataRCALz =   -153.88 )
      parameter ( mcRCALz   =   -152.13 )
      real datahesz,mchesz, hesz
      parameter ( datahesz  =   -155.21 ) ! best position for no bias
      parameter ( mchesz    =   -154.14 ) ! mc hes Z
      real dataFhesz,mcFhesz, fhesz
      parameter ( dataFhesz  =   227.2 )
      parameter ( mcFhesz    =   227.2 )
      real datasrtdz, mcsrtdz, srtdz
      parameter ( datasrtdz =  -149.31 )  ! same as HES Z shift
      parameter ( mcsrtdz   =  -148.25 )  ! midpoint of 2 srtd planes

      integer whatToUse
      integer useCAL, useHES, useSRTD, useVCTRHL, useVCTPAR
      parameter ( useCAL    = 1 )
      parameter ( useHES    = 2 )
      parameter ( useSRTD   = 3 )
      parameter ( useVCTRHL = 4 )
      parameter ( useVCTPAR = 5 )

      real calradius
      real modpos
      real hesgap
C     leave a gap of 2cm on BOTH SIDES of Module crack
      parameter ( hesgap = 2.0 )   
      real SRTDedge
C     Use SRTD away from the edge by SRTDedge cm
      parameter ( SRTDedge = 2.0 )   

      real MAXMIP
      PARAMETER ( MAXMIP = 1000. ) ! stay away from high default values. 


      if ( mode.eq.1.and.cal(1) .eq. savecal(1) .and. cal(2) .eq.
     &     savecal(2) ) then 
C
C  Don't re-calculate if the user just asks again for best the best theta
C
         best = lastBestTheta
         return
      endif
      
C
C  Save cal position as a record of the event
C
      savecal(1) = cal(1)
      savecal(2) = cal(2)

      if (montecarlo) then
         rcalz = mcrcalz
         hesz  = mchesz
         fhesz = mcfhesz
         srtdz  = mcsrtdz
      else
         rcalz = datarcalz
         hesz  = datahesz
         fhesz = datafhesz
         srtdz  = datasrtdz
      endif

      WhatToUse = useCAL  ! start by using cal
      hesexpected = 0
      hesfound  = 0
      calradius = sqrt(cal(1)**2+cal(2)**2)

C -------------------------------
C  See if we can use HES
C -------------------------------

C
C  IN RCAL use HES over CAL when not in module cracks and
C  the position is more than 20 cm from RCAL beampipe
C  (near the beampipe the hes diodes end and the SRTD
C  takes over anyway.
C
      if ( cal(3).lt.-140 .and. calRadius.gt.20.) then
C     in the RCAL more than 20 cm from beampipe.

         modpos = abs(cal(1)) - 20.33/2. + hesgap
         if ( modpos - int(modpos/20.33)*20.33 .gt. 2*hesgap .or.
     &        modpos.lt.0 ) then
C
C     The cal says that we are NOT IN THE MODULE CRACK so I expect there
C     to be a HES MIP.
C
            hesexpected = 1
            if (hes(4).ge.5..and.hes(4).lt.MAXMIP) then
               hesfound = 1
               WhatToUse = UseHES
            endif
         endif ! not in crack
      endif ! in rcal
C
C  In FCAL use the hes whenever available
C
      if ( cal(3).gt.220.and.hes(4).gt.5.and.hes(4).lt.MAXMIP) then
         WhatToUse = UseHES
      endif

C -------------------------------
C  See if we can use SRTD
C
C  require an cal position to be in RCAL and the SRTD mips to be o.k.
C  then require that the electron is within SRTDedge cm of the srtd edge.
C
C -------------------------------
      if ( cal(3).lt.-140 .and.srtd(4).gt.0.4.and.srtd(4).lt.MAXMIP )
     &     then
         if (srtd(1).gt.(-34+srtdedge).and.srtd(1).lt.(34-srtdedge).and
     &        .(srtd(1).gt.(-10-srtdedge).and.srtd(2).gt.(-28+srtdedge)
     &        .or.srtd(1).le.(-10-srtdedge).and.srtd(2).gt.-40+srtdedge)
     &        .and.(srtd(1).ge.(10+srtdedge).and.srtd(2).lt.(40-srtdedge
     &        ).or.srtd(1).lt.(10+srtdedge).and.srtd(2).lt.28-srtdedge)
     &        )then
            WhatToUse = UseSRTD
         endif
      endif

C -------------------------------
C  See if we can use VCTRHL dip angle
C
C  We must be more than 80 cm away from 0,0 for there to be
C  enough stereo hits for the track to give good resolution.
C -------------------------------
      if ( calRadius.gt.80) then
         WhatToUse = UseVCTRHL
      endif

C -------------------------------
C  See if we can use VCTPAR  angle
C
C  VCTPAR angle must be in range and the VCTRHL dip requirements
C  must already be satisfied.
C -------------------------------
c      if ( tvctpar.gt.0.and.tvctpar.lt.3.2.and.useVCTRHL) then
c         WhatToUse = UseVCTPAR
c      endif

C -------------------------------------
C  Now I know what to use.
C  calculate the best angle.
C
C  if mode = 2 then calculate all angles.
C -------------------------------------


C
C  Start with Elecpo
C
      if ( WhatToUse .eq. UseCal .or. mode.eq.2) then
         if ( cal(3).lt.-140 ) then
C     ...  USE SHIFTED RCAL
            thetacal = atan2(sqrt(cal(1)**2+cal(2)**2),(rcalz-vtx(3)))
         else
            thetacal = atan2(sqrt(cal(1)**2+cal(2)**2),(cal(3)-vtx(3)))
         endif
      endif

C
C  Calculate HES theta
C
      if ( WhatToUse .eq. UseHES .or. (mode.eq.2.and.hes(4).ge.5.and
     $     .hes(4).lt.MAXMIP)) then
         if ( cal(3).lt.-140 ) then
C     ...  USE SHIFTED RCAL
            thetahes = atan2(sqrt(hes(1)**2+hes(2)**2),(hesz-vtx(3)))
         else
            thetahes = atan2(sqrt(hes(1)**2+hes(2)**2),(fhesz-vtx(3)))
         endif
      endif

C
C  Calculate SRTD Theta
C
      if ( WhatToUse .eq. UseSRTD .or. (mode.eq.2.and.srtd(4).gt.0.4.and
     &     .srtd(4).lt.MAXMIP)) then
            thetasrtd = atan2(sqrt(srtd(1)**2+srtd(2)**2),
     &        (srtdz-vtx(3)))
      endif


C
C  Calculate VCTRHL Theta
C

      if ( WhatToUse .eq. UseVCTRHL .or.(mode.eq.2.and.calradius.gt.50))
     &     then
         thetavctrhl = tdip
         if (thetaVCTRHL.lt.0) thetaVCTRHL = 3.14159 + thetaVCTRHL
      endif

C
C  Calculate VCTPAR Theta
C  (not needed since VCTPAR is already an angle)



      best = -1.

      if (whatToUse .eq. useCal) best = thetacal
      if (whatToUse .eq. useHES) best = thetaHES
      if (whatToUse .eq. useSRTD) best = thetaSRTD
      if (whatToUse .eq. useVCTRHL) best = thetaVCTRHL
      if (whatToUse .eq. useVCTPAR) best = tvctpar

      lastbestTheta = best

      return
      end
      
