C     ===================================================
      SUBROUTINE GETDA (EE,EP,YJB,Q2JB,XDA,Q2DA,YDA,TH_E)
C     ===================================================

      IMPLICIT NONE

      REAL EE,EP,YJB,Q2JB,XDA,Q2DA,YDA,TH_E,A,B,COSG,G,SING,X0,
     &     COSTH,SINTH,BETA,DELTA

      REAL NO_USE
      DATA NO_USE / -999.99 /

      A    = Q2JB*(1.0-YJB)
      B    = 4.0*(EE**2)*(YJB**2)

      IF ( (A+B).EQ.0 ) THEN
        YDA  = NO_USE
        XDA  = NO_USE
        Q2DA = NO_USE
        GOTO 100
      ENDIF

      COSG = (A-B)/(A+B)
      G    = ACOS(COSG)
      SING = SIN(G)

      X0     = EE/EP
      COSTH  = COS(TH_E)
      SINTH  = SIN(TH_E)
      BETA   = SING + SINTH - SIN(TH_E+G)
      DELTA  = SING + SINTH + SIN(TH_E+G)

      IF (BETA.NE.0) THEN
        YDA    = ( SINTH*(1-COSG) )/ BETA
        Q2DA   = ( 4*(EE**2)*SING*(1+COSTH) )/ BETA
        XDA    = (X0*DELTA) / BETA
      ELSE
        YDA  = NO_USE
        XDA  = NO_USE
        Q2DA = NO_USE
      ENDIF

100   CONTINUE

      RETURN
      END


