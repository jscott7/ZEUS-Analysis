C     ------
C  Parameter for array MAX
      INTEGER LIMIT_O1CNVD
            PARAMETER( LIMIT_O1CNVD = 999 )
C
C  O1CNVD      : GFLT component FLT data; unit is according to
C                offline-convention
C  MAX_O1CNVD  : all O1CNVD(IPNT) = NULLDT for IPNT .GT. MAX_O1CNVD
C  TOT_O1CNVD  : number of available( .NE. NULLDT) component data
C
      INTEGER MAX_O1CNVD, TOT_O1CNVD, O1CNVD
      COMMON /O1CNVD/ MAX_O1CNVD,
     +                TOT_O1CNVD,
     +                O1CNVD(LIMIT_O1CNVD)
C
