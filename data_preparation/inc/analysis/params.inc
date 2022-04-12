*-->  initial energy of proton and e+ and other stuff
      REAL       Ep0         , eE0
      PARAMETER( Ep0=820.001 , eE0=27.52 )

      Real       M2p
      Parameter( M2p = 0.9382723**2 )

      Real       twoeE0
      PARAMETER( twoeE0 = 2.0*eE0 )

      Real       fourEp0eE0
      Parameter( fourEp0eE0 = 4.0*Ep0*Ee0 )

      Real       x0
      PARAMETER( x0 = eE0/Ep0 )

      Real       MIN_Q2        , MAX_Q2
      Parameter( MIN_Q2 = 1.0  , MAX_Q2 = 1.E5  ) 

      Real       MIN_Y         , MAX_Y 
      Parameter( MIN_Y = 1.E-4 , MAX_Y = 0.9999 )

      Real       MIN_X         , MAX_X
      Parameter( MIN_X = 1.E-5 , MAX_X = 0.9999 )

*-->  some usefull variables
      Integer    I,K
      Real       pte2

*-->  geometric variables
      Real       Cgamma1, Cgamma2 , Cgamma3, Cgamma4
      REAL       Ctheta, Cgmtru
      REAL       Pi, TwoPi, PiBy2, PiBy4
      PARAMETER (Pi = 3.14159265358979323846)
      PARAMETER (TwoPi = 2.*Pi)
      PARAMETER (PiBy2 = Pi/2.0) 
      PARAMETER (PiBy4 = Pi/4.0) 

      Real       piMass
      Parameter (piMass = 0.13956995)
      Real       xpi2
      Parameter (xpi2   = piMass*piMass)

