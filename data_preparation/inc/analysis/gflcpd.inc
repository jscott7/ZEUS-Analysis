C     ------
C Common to store component FLT data
C----
C  Parameter for array MAX
      INTEGER LIMIT_O1CPDT
            PARAMETER( LIMIT_O1CPDT = 999 )
C
C  O1CPDT      : GFLT component FLT data based upon O1RCDT_Rcdat
C  MAX_O1CPDT  : all O1CPDT(IPNT) = NULLDT for IPNT .GT. MAX_O1CPDT
C  TOT_O1CPDT  : number of available( .NE. NULLDT) component data
C
      INTEGER MAX_O1CPDT, TOT_O1CPDT, O1CPDT
      COMMON /O1CPDT/ MAX_O1CPDT,
     +                TOT_O1CPDT,
     +                O1CPDT(LIMIT_O1CPDT)
C
C     ------
