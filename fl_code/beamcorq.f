      REAL FUNCTION BEAMCORQ( Gamma_PT , Q2_PT, Icheck, Ierr)
      Implicit None
C     ---------------------------------------------------------------
C     Input: Gamma_PT : gamma, calculated with the PT method ...
C                      ... in RADIANS !!
C            Q2_PT    : Let's make the correction Q2 dependent...
C            Icheck   2: Default  !
C                     1: Lower Correction
C                     3: Higher Correction
C     Output: Ierr
C             1: Gamma not in Radians??
C             2: Q2 out of range
C
C     Purpose: Correct MonteCarlo for wrong beampipe simulation in 96/97
C     Author: N. Tuning 08/04/00
C
C     Changed from beamcor.f to beamcorq.f 15/04/00
C     Version 2.2 16/04/00
C     --------------------
C
C     How does the correction factor look?? Plot this:
C     PAW> opt logx
C     PAW> null 0.01 3.2 0.5 1.5
C     PAW> fun/pl beamcorq.f(x,20.,2,0) 0.01 3.2 s
C     PAW> atit 'Gamma?PT! (rad)' 'Beampipe Corr. Factor'
C     ---------------------------------------------------------------
      Real    BEAMCOR1, BEAMCOR2, BEAMCOR3
      Real    Gamma_PT, Q2_PT, X
      Integer Icheck, Ierr

      Real Shift_Lowg, Scale_Lowg, Slope_Lowg
      Real Shift_Higg, Scale_Higg, Slope_Higg
      Data Shift_Lowg, Scale_Lowg, Slope_Lowg / 0.85, 14., 4.5 /
      Data Shift_Higg, Scale_Higg, Slope_HIgg / 0.25,  4., 5.  /
      Real Q2shift
      Real Shift_Lowg2, Shift_Higg2
      Logical First
      Data    First/.TRUE./
C     --------------------------------------------------------------
      If (First) Then
         First=.false.
         Write(*,*) '=========================================='
         Write(*,*) ' This is Beampipe Correction version v2.2 '
         Write(*,*) '                                          '
         Write(*,*) ' Date 16/04/2000 (N.Tuning)               '
         Write(*,*) ' Only to be used for 1996/1997 analysis   '
         Write(*,*) '                                          '
         Write(*,*) '=========================================='
      EndIf
      

      Ierr = 0
      BEAMCORQ = 1.
C     =============
      IF (Gamma_PT.LE.0.) THEN
         Ierr = 1
         Write(*,*) 'Gamma_Pt le 0 !! (beamcorq.f):',Gamma_PT
         RETURN
      ELSEIF (Gamma_PT.GT.3.2) THEN
         Ierr = 1
         Write(*,*) 'Gamma_Pt > Pi !! (beamcorq.f):',Gamma_PT
         Write(*,*) '(gamma should bne in radians, not in degrees...)'
         RETURN
      ENDIF

      IF (Q2_PT.LE.0.) THEN
         Ierr = 2
         Write(*,*) 'Q2 le 0. !!(beamcorq.f):',Q2_PT
         RETURN
      ENDIF

C     ===================
      X = log10(Gamma_PT)
C     ===================
      Q2shift = 0.002*(log10(Q2_PT)-1.6)

      Shift_Higg2 = Shift_Higg+Q2shift
      Shift_Lowg2 = Shift_Lowg+Q2shift

      BeamCor1 = 
     +     Scale_Higg*( 
     +     ( exp( Slope_Higg*(x+Shift_Higg2))
     +     - exp(-Slope_Higg*(x+Shift_Higg2)))    /
     +     ( exp( Slope_Higg*(x+Shift_Higg2))
     +     + exp(-Slope_Higg*(x+Shift_Higg2)))  -1 )
     +
     +     -1.*Scale_Lowg*( 
     +     ( exp( Slope_Lowg*(x+Shift_Lowg2))
     +     - exp(-Slope_Lowg*(x+Shift_Lowg2)))      /
     +     ( exp( Slope_Lowg*(x+Shift_Lowg2))
     +     + exp(-Slope_Lowg*(x+Shift_Lowg2)))    -1 )

      BeamCor2 = 1. + (1./100)*BeamCor1

C---
C Systematic Check:
C---
      If (Icheck.eq.1) Then
         BeamCor3 = Beamcor2 + 0.01/(Gamma_PT+0.09)
      ElseIf (Icheck.eq.3) Then
         BeamCor3 = Beamcor2 - 0.01/(Gamma_PT+0.09)
      Else
         BeamCor3 = Beamcor2
      EndIf

      If (BeamCor3.gt.0.) Then
C     ==========================
         BeamCorQ = 1./BeamCor3
C     ===========================
      Else
         Write(*,100) 
     +        Beamcor1,Beamcor2,BeamCor3,
     +        gamma_pt,Q2_PT,Q2shift,
     +        Icheck
      EndIf
 100  FORMAT('ERROR beamcorq.f beamcor3:',6F12.5,I3)
      End



