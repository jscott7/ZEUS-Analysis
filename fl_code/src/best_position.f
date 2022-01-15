      subroutine best_position(besttheta,bestx,besty,bestz,cal,srtd,hes
     &     ,Tend,vtx,montecarlo,run,whatToUse,Isys)
C
C best_position returns the best reconstructed position and angle for
C     the scattered electron given the CAL, SRTD, HES and track endpoint
C     for the electron and the measured vertex.
C
C  For Rcal electrons (Z<-140) the Z is projected onto the CAL plane.
C
C  This routine expects that the modified CCCXYZ was run in the 
C  eaze job to shift RCAL Z among other shifts.
C
C OUTPUT: 
C  real besttheta = the resulting best angle for the electron
C  real bestx,besty,bestz = the resulting best position projected onto 
C                           an Rcal,Fcal plane or Bcal cylinder
C
C  real srtd(1)   = the shifted SRTD x position
C
C
C INPUT:
C    Real cal(4)  
C         cal(1) = x cal from elecpo
C         cal(2) = y cal from elecpo
C         cal(3) = z cal from elecpo ( not used)
C         cal(4) = electron energy (not used by best_theta)
C
C    Real srtd(4)  
C         srtd(1) = x from srtdelec (shifted after call)
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
C    Real Tend(6)  
C         tend(1) = x from Track
C         tend(2) = y from Track
C         tend(3) = z from track (not used)
C         tend(4) = Track Momentum 
C                   (if Momentum > 5 then Tend can be used)
C         tend(5) = track DCA
C                   (DCA must be < 10 for this to be used)
C         tend(6) = track dip value from VCTRHL
C                   ( tdip < 100 for this to be used )
C                   ( If you do not want to use tracking angle set this
C                      value to 1000 )
C
C    Real vtx(3) = vertex (x,y,z) position
C
C    Logical Montecarlo = True if montecarlo
C
C    Integer run  = Data run number
C
C Author: Mike Wodarczyk  13-June-1998   Version 1.1
C
C Changes: Jonathan Scott 8-February-2000
C Include code to catch MC ntuple errors.
C
      Implicit None

      Integer Isys

      Real cal(4),srtd(4),hes(4),tend(6),vtx(3)
      Logical Montecarlo
      Real bestx, besty,bestz,besttheta
      integer run

      real RCAL_FACE
      parameter ( RCAL_FACE   =   -152.13 )
      real BCAL_FACE
      parameter ( BCAL_FACE   =   125.669 )
      real FCAL_FACE
      parameter ( FCAL_FACE   =   226.13  )

c added 22/06/99 for check_no

C_AB      include '../inc/xq2bins_sys.inc'
      real    srtd_sh

      real dataRCALz,mcRCALz, rcalz
C      parameter ( dataRCALz =   -153.88 )
      parameter ( dataRCALz =   -152.13 )
C      parameter ( dataRCALz =   -153.03 )
      parameter ( mcRCALz   =   -152.13 )
      real datahesz,mchesz, hesz
      parameter ( datahesz  =   -155.04 ) ! best position for no bias
      parameter ( mchesz    =   -154.14 ) ! mc hes Z
      real dataFhesz,mcFhesz, fhesz
      parameter ( dataFhesz  =   227.2 )
      parameter ( mcFhesz    =   227.2 )
      real datasrtdz, mcsrtdz, srtdz
      parameter ( datasrtdz =  -149.15 )  ! same as HES Z shift
      parameter ( mcsrtdz   =  -148.25 )  ! midpoint of 2 srtd planes

      integer whatToUse
      integer useCAL, useHES, useSRTD, useTrack
      parameter ( useCAL    = 1 )
      parameter ( useHES    = 2 )
      parameter ( useSRTD   = 3 )
      parameter ( useTrack = 4 )

      real calradius, radius, Caltheta
      real modpos
      real hesgap
C     leave a gap of 2cm on BOTH SIDES of Module crack
      parameter ( hesgap = 2.0 )   
      real SRTDedge
C     Use SRTD away from the edge by SRTDedge cm
      parameter ( SRTDedge = 2.0 )   

      real MAXMIP
      PARAMETER ( MAXMIP = 1000. ) ! stay away from high default values. 
      real MAXTDIP
      PARAMETER ( MAXTDIP = 100. ) ! stay away from high default values. 

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
      if (cal(3).lt.-140.) then 
	  if ((cal(1)**2+cal(2)**2).ge.0.) then
          calradius = sqrt(cal(1)**2+cal(2)**2)
	  else
	    write(*,*) 'ERROR: NTUPLE INPUT BAD'
	    return
	  end if
      else
       calradius = 999.
      endif


      Caltheta = atan2(sqrt(cal(1)**2+cal(2)**2),(cal(3)-vtx(3)))


C -------------------------------
C  See if we can use HES
C -------------------------------

C
C  IN RCAL use HES over CAL when not in module cracks and
C  the position is more than 20 cm from RCAL beampipe
C  (near the beampipe the hes diodes end and the SRTD
C  takes over anyway.
C
C  Exclude HES in Rcal module 12 since the position is ugly.
C
 
      if ( cal(3).lt.-140 .and. calRadius.gt.20..and. abs(cal(1)).gt.10.
     &     ) then

C     in the RCAL more than 20 cm from beampipe.

         modpos = abs(cal(1)) - 20.33/2. + hesgap
         if ( modpos - int(modpos/20.33)*20.33 .gt. 2*hesgap .or.
     &        modpos.lt.0 ) then
C
C     The cal says that we are NOT IN THE MODULE CRACK so I expect there
C     to be a HES MIP.
C
            if (hes(4).ge.5..and.hes(4).lt.MAXMIP) then
               WhatToUse = UseHES
            endif
         endif ! not in crack
      endif ! in rcal
C
C  In FCAL use the hes whenever available (taken out 02/06/99)
C
c$$$      if ( cal(3).gt.220.and.hes(4).gt.5.and.hes(4).lt.MAXMIP) then
c$$$         WhatToUse = UseHES
c$$$      endif

C -------------------------------
C  See if we can use SRTD
C
C  require an cal position to be in RCAL and the SRTD mips to be o.k.
C  then require that the electron is within SRTDedge cm of the srtd edge.
C
C -------------------------------
      if ( cal(3).lt.-140 .and.srtd(4).gt.0.4.and.srtd(4).lt.MAXMIP )
     $     then
         
         if (srtd(1).gt.(-34+srtdedge).and.srtd(1).lt.(34-srtdedge).and
     $        .(srtd(1).gt.(-10-srtdedge).and.srtd(2).gt.(-28+srtdedge)
     $        .or.srtd(1).le.(-10-srtdedge).and.srtd(2).gt.-40+srtdedge)
     $        .and.(srtd(1).ge.(10+srtdedge).and.srtd(2).lt.(40-srtdedge
     $        ).or.srtd(1).lt.(10+srtdedge).and.srtd(2).lt.28-srtdedge)
     $        )then
            WhatToUse = UseSRTD
         endif

         if (.not.montecarlo.and. run.ge.25296 ) then
            srtd(1) = srtd(1) + 0.2 ! shift X
         endif
         if ( .not.montecarlo ) then
C     
C     Data crack wider by 0.1 cm  on left and right
C     
c taken out 25/06/99

c            if ( (srtd(1).gt.10.and.srtd(2).gt.0 ).or. 
c     &           (srtd(1).gt.-10.and.srtd(2).lt.0)) then
c               srtd(1) = srtd(1) + 0.05 ! shift X 
c            else
c               srtd(1) = srtd(1) - 0.05 ! shift X 
c            endif

c
c 25/06/99: shift  SRTD quadrants seperately
c
c         II                  I
c       +4 mm               0.0 mm
c
c         III                IV
c       +2.5 mm             -1.5 mm
c

c
c 08/09/99: shifts slightly changed to numbers cited in thesis
c
c         II                  I
c       +4 mm               -0.5 mm
c
c         III                IV
c       +2.5 mm             -2.5 mm
c




         if (srtd(1).gt.10.0.and.srtd(2).gt.-4.3)  srtd_sh = -0.05    
         if (srtd(1).lt.10.0.and.srtd(2).gt.3.7) srtd_sh = 0.4
         if (srtd(1).lt.-10.0.and.srtd(2).lt.3.7)  srtd_sh = 0.25    
         if (srtd(1).gt.-10.0.and.srtd(2).lt.-4.3) srtd_sh = -0.25

c --------------------------------------------------
c check 53: move SRTD by additional 2 mm in y (data)
c HAS TO BE DONE AFTER CHECK OF FIDUCIAL VOLUME!!!
c --------------------------------------------------

c          if (check_no.eq.53) srtd_sh = 0.7          
          if (Isys.eq.53) srtd_sh = srtd_sh + 0.05


c ---------------------------------------------  
c check 54: move SRTD by 2 mm less in y (data)
c ---------------------------------------------

c          if (check_no.eq.54) srtd_sh = 0.3
          if (Isys.eq.54) srtd_sh = srtd_sh - 0.05

c$$$            if  ( srtd(1).lt.-10.or.(srtd(1).lt.10.and.
c$$$     &            srtd(2).gt.0)) then 
c$$$               srtd(2) = srtd(2) + srtd_sh
c$$$            endif

          srtd(2) = srtd(2) + srtd_sh
  
         endif

      endif

C -------------------------------
C  See if we can use Track Endpoint 
C
C  We must be more than 80 cm away from 0,0 for there to be
C  enough stereo hits for the track to give good resolution.
C -------------------------------
      if ( calRadius.gt.80.0 .and. caltheta.gt.0.3) then
         WhatToUse = UseTrack
      endif

C
C     Start with Elecpo
C
      besttheta = -1.

      if ( WhatToUse .eq. UseCal ) then
            bestx = cal(1)
            besty = cal(2)
            bestz = cal(3)
            besttheta = atan2(sqrt(cal(1)**2+cal(2)**2),(cal(3)-vtx(3)))
      endif

C
C  Calculate HES theta
C
      if ( WhatToUse .eq. UseHES ) then
         if ( cal(3).lt.-140 ) then
C     ...  USE SHIFTED RCAL
            bestx = hes(1)
            besty = hes(2)
            bestz = hesz
            besttheta = atan2(sqrt(hes(1)**2+hes(2)**2),(hesz-vtx(3)))
         else
            bestx = hes(1)
            besty = hes(2)
            bestz = fhesz
            besttheta = atan2(sqrt(hes(1)**2+hes(2)**2),(fhesz-vtx(3)))
         endif
      endif

C
C  Calculate SRTD Theta
C
      if ( WhatToUse .eq. UseSRTD ) then
         bestx = srtd(1)
         besty = srtd(2)
         bestz = srtdz
         besttheta = atan2(sqrt(srtd(1)**2+srtd(2)**2),
     $        (srtdz-vtx(3)))
      endif

C
C  Calculate VCTRHL Theta
C
      if ( WhatToUse .eq. UseTrack) then
         bestx = Tend(1)
         besty = Tend(2)
         bestz = Tend(3)
          if (besttheta.lt.0) besttheta = 3.14159 + besttheta
      endif

C
C  Now project all x,y,z onto planes at the RCAL FACE or Bcal cylinder
c     or FCAL face
C

      if ( bestz .lt. -140 ) then
C RCAL
         bestx = bestx * (RCAL_FACE - vtx(3)) / (bestz - vtx(3))
         besty = besty * (RCAL_FACE - vtx(3)) / (bestz - vtx(3))
         bestz = RCAL_FACE
      else if (bestz.lt.220 ) then
C BCAL
         radius = sqrt(bestx**2+besty**2) 
         bestx = bestx*BCAL_FACE / radius
         besty = besty*BCAL_FACE / radius
         bestz = vtx(3) + (bestz-vtx(3))*BCAL_FACE / radius
      else
C FCAL
         bestx = bestx * (FCAL_FACE - vtx(3)) / (bestz - vtx(3))
         besty = besty * (FCAL_FACE - vtx(3)) / (bestz - vtx(3))
         bestz = FCAL_FACE
      endif
      return
      end
      
