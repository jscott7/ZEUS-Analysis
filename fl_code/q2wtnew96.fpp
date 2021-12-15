C
C     ===========================================================
      SUBROUTINE Q2WTNEW96(Q2,WTQ2,IN)
C     ===========================================================
C
      Implicit None

#include "local.inc"

      LOGICAL Q2WTFIRST
      DATA  Q2WTFIRST /.TRUE./

      REAL WTQ2,WTDQ2,WTNDQ2,Q2, lumc1, lumc2, lumc3, lumc4
      INTEGER IN

      IF (Q2WTFIRST.EQV..TRUE.) THEN
       LUMC1=LUMCND1+LUMCD1
       LUMC2=LUMCND1+LUMCD1+LUMCND2+LUMCD2
       LUMC3=LUMCND1+LUMCD1+LUMCND2+LUMCD2+LUMCND3+LUMCD3
       LUMC4=LUMCND1+LUMCD1+LUMCND2+LUMCD2+LUMCND3+LUMCD3
     & +LUMCND4+LUMCD4
       WRITE(6,*) ' *******************************'       
       WRITE(6,*) ' '       
       WRITE(6,*) ' LUMI DATA :',LUMDA
       WRITE(6,*) ' LUMI MC :'
       WRITE(6,*) LUMC1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(6,*) LUMC2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(6,*) LUMC3,'for',Q2MC3,'< Q2 <',Q2MC4
       WRITE(6,*) LUMC4,'for',Q2MC4,'< Q2'
       WRITE(6,*) ' '       
       WRITE(6,*) ' '       
       WRITE(6,*) ' LUMI FOR NONDIF MC :'
       WRITE(6,*) LUMCND1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(6,*) LUMCND2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(6,*) LUMCND3,'for',Q2MC3,'< Q2 <',Q2MC4
       WRITE(6,*) LUMCND4,'for',Q2MC4,'< Q2'
       WRITE(6,*) ' '       
       WRITE(6,*) '    LUMI FOR DIF MC :'
       WRITE(6,*) LUMCD1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(6,*) LUMCD2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(6,*) LUMCD3,'for',Q2MC3,'< Q2 <',Q2MC4
       WRITE(6,*) LUMCD4,'for',Q2MC4,'< Q2'
       WRITE(6,*) ' '       
       WRITE(6,*) ' *******************************'       
       WRITE(999,*) ' *******************************'       
       WRITE(999,*) ' '       
       WRITE(999,*) ' LUMI DATA :',LUMDA
       WRITE(999,*) ' LUMI MC :'
       WRITE(999,*) LUMC1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(999,*) LUMC2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(999,*) LUMC3,'for',Q2MC3,'< Q2 <',Q2MC4
       WRITE(999,*) LUMC4,'for',Q2MC4,'< Q2'
       WRITE(999,*) ' '       
       WRITE(999,*) ' '       
       WRITE(999,*) ' LUMI FOR NONDIF MC :'
       WRITE(999,*) LUMCND1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(999,*) LUMCND2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(999,*) LUMCND3,'for',Q2MC3,'< Q2 <',Q2MC4
       WRITE(999,*) LUMCND4,'for',Q2MC4,'< Q2'
       WRITE(999,*) ' '       
       WRITE(999,*) '    LUMI FOR DIF MC :'
       WRITE(999,*) LUMCD1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(999,*) LUMCD2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(999,*) LUMCD3,'for',Q2MC3,'< Q2 <',Q2MC4
       WRITE(999,*) LUMCD4,'for',Q2MC4,'< Q2'
       WRITE(999,*) ' '       
       WRITE(999,*) ' *******************************'       

       Q2WTFIRST=.FALSE.
      ENDIF
    
      WTQ2=0.0 

      IF (Q2.GT.Q2MC1.AND.Q2.LE.Q2MC2) THEN
       WTDQ2=LUMDA/LUMCD1
       WTNDQ2=LUMDA/LUMCND1

      ELSEIF (Q2.GT.Q2MC2.AND.Q2.LE.Q2MC3) THEN
       WTDQ2=LUMDA/(LUMCD1+LUMCD2)
       WTNDQ2=LUMDA/(LUMCND1+LUMCND2)

      ELSEIF (Q2.GT.Q2MC3.AND.Q2.LE.Q2MC4) THEN
       WTDQ2=LUMDA/(LUMCD1+LUMCD2+LUMCD3)
       WTNDQ2=LUMDA/(LUMCND1+LUMCND2+LUMCND3)

      ELSEIF (Q2.GT.Q2MC4) THEN
       WTDQ2=LUMDA/(LUMCD1+LUMCD2+LUMCD3+LUMCD4)
       WTNDQ2=LUMDA/(LUMCND1+LUMCND2+LUMCND3+LUMCND4)
      ENDIF
      
      IF (IN.GE.1.AND.IN.LE.NNODIFFILES) THEN
       WTQ2=WTNDQ2*(1.0-FRACDIF)

      ELSEIF (IN.GT.NNODIFFILES.AND.IN.LE.NMCFILES) THEN
       WTQ2=WTDQ2*FRACDIF

      ELSEIF (IN.LT.1.OR.IN.GT.NMCFILES) THEN
       WRITE (6,*) ' '
       WRITE (6,*) '  ERROR IN WTQ2  '
       WRITE (6,*) ' '
       WRITE (999,*) ' '
       WRITE (999,*) '  ERROR IN WTQ2  '
       WRITE (999,*) ' '

      ENDIF

      RETURN
      END

C                                                                           


