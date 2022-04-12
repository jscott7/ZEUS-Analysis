C Common for Component FLT data made with ZGANA
C ( GFLCPD in o1recon is the corresponding sequence. )
C
C  Parameter for array MAX
      INTEGER LIMIT_O1RCEX
            PARAMETER( LIMIT_O1RCEX = 999 )
C
C  O1RCEX      : GFLT component FLT data made with ZGANA from readout data
C  MAX_O1RCEX  : all O1RCEX(IPNT) = NULLDT for IPNT .GT. MAX_O1RCEX
C  TOT_O1RCEX  : number of available( .NE. NULLDT) component data
C
      INTEGER MAX_O1RCEX, TOT_O1RCEX, O1RCEX
      COMMON /O1RCEX/ MAX_O1RCEX,
     +                TOT_O1RCEX,
     +                O1RCEX(LIMIT_O1RCEX)
C
