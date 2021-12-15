      subroutine dead_material_energy(EnergyOut,EnergyIn,cal,srtd,pres
     $     ,montecarlo,mode,scalefactor)

C
C  Gives the presampler or srtd corrected energy for the electron
C
C  Input:
C     real EnergyIn = Energy of electron before presampler or
C                     SRTD corrections. This may be the same or different from
C                     cal(4) depending on if other (non-uniformity) corrections 
C                     have been applied. 
C
C     real cal(4):
C       cal(1) = elecpo X of electron (not used)
C       cal(2) = elecpo Y of electron (not used)
C       cal(3) = elecpo Z of electron 
C       cal(4) = Cal energy (not used)
C
C     real srtd(4):
C      srtd(1) = srtd X 
C      srtd(2) = srtd Y 
C      srtd(3) = srtd Z (not used)
C      srtd(4) = mips from SRTDELEC (sum of x and y planes)
C
C     real pres = presampler mips sum from  3x3 tiles. 
C         This is the second element of the presampler array from prclus.
C
C     integer mode = way to perform correction for systematic checks
C         mode = 0  Apply nominal correction
C         mode = 1  Shift data correction up by systematic error and
C                         mc down by systematic error.
C         mode = 2  Shift data correction down by systematic error and
C                         mc up by systematic error.
C
C     Logical Montecarlo = True if montecarlo. False if data.
C
C  Output:
C     real EnergyOut = Electron energy after presampler or SRTD corrections. 
C
C
C     Version 1.1  July 7, 1998 Mike Wodarczyk
C         Changed the presampler and SRTD correction values. (code is the same).
C         The new values are based on KP and compton corrections fitting the 
C         peak of (Cal+corr*mips) / Predicted and choosing corr so the peak=1.
C
      implicit none
	real scalefactor
      real EnergyOut,EnergyIn
      real cal(4),srtd(4), pres
      logical montecarlo
      integer mode

      integer whatToUse
      integer useCAL, useSRTD, usePRES
      parameter ( useCAL    = 1 )
      parameter ( useSRTD   = 2 )
      parameter ( usePRES   = 3 )

      real MAXMIP
      PARAMETER ( MAXMIP = 1000. ) ! stay away from high default values. 

      real McSrtdCorr
      parameter ( McSrtdCorr = 0.039 ) 
      real DataSrtdCorr
      parameter ( DataSrtdCorr = 0.039 ) 
      real McSrtdCorrErr
      parameter ( McSrtdCorrErr = 0.005 ) 
      real DataSrtdCorrErr
      parameter ( DataSrtdCorrErr = 0.005 ) 

      real McPresCorr
      parameter ( McPresCorr = 0.0790 ) 
      real DataPresCorr
      parameter ( DataPresCorr = 0.0654 ) 
      real McPresCorrErr
      parameter ( McPresCorrErr = 0.005 ) 
      real DataPresCorrErr
      parameter ( DataPresCorrErr = 0.005 ) 

      real SRTDedge
C     Use SRTD away from the edge by SRTDedge cm
      parameter ( SRTDedge = 2.0 )   


      real PresCorr,SRTDCorr, PresErr, SRTDErr

      if ( montecarlo) then
         PresCorr = McPresCorr
         SRTDCorr = McSRTDCorr
C  for MC systematic, shift the correction the opposite way relative to
C     the data shift.
         PresErr  = -1.*McPresCorrErr  
         SRTDErr  = -1.*McSRTDCorrErr
      else
         PresCorr = DataPresCorr
         SRTDCorr = DataSRTDCorr
         PresErr  = DataPresCorrErr
         SRTDErr  = DataSRTDCorrErr
      endif

      if ( mode.eq. 1) then
         PresCorr = PresCorr + PresErr
         SRTDCorr = SRTDCorr + SRTDErr
      else if ( mode .eq. 2) then
         PresCorr = PresCorr - PresErr
         SRTDCorr = SRTDCorr - SRTDErr
      else if (mode .ne. 0) then
         write(*,*) 'Error in dead_material_energy.'
         write(*,*) ' Mode = ',WhatToUse
         stop
      endif


      WhatToUse = UseCAL


C -------------------------------
C  See if we can use Presampler
C 
C  Require an RCAL electron that it has mips in valid range
C  I saw no events between 0 and 0.025 mips in data.  I use
C  0.025 as a lower cut off so there are no "0" rounding problems.
C -------------------------------
      if ( cal(3).lt.-140 .and. pres.gt.0.025 .and. pres.lt.MAXMIP) then
         WhatToUse  = UsePRES
      endif


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
      endif

      if (WhatToUse .eq. useCal) then
         EnergyOut = EnergyIn         ! by default use the cal
      else if (WhatToUse .eq. usePres) then
         if (montecarlo) then
            EnergyOut = EnergyIn + pres*McPresCorr*scalefactor
         else
            EnergyOut = EnergyIn + pres*DataPresCorr*scalefactor
         endif
      else if (WhatToUse .eq. useSRTD) then
         if (montecarlo) then
            EnergyOut = EnergyIn + srtd(4)*McSrtdCorr*scalefactor
         else
            EnergyOut = EnergyIn + srtd(4)*DataSrtdCorr*scalefactor
         endif
      else
         write(*,*) 'Error in dead_material_energy.  coding error.'
         write(*,*) ' WhatToUse = ',WhatToUse
         stop
      endif

      return
      end


