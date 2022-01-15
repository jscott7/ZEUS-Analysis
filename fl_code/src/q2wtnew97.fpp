C
C     ===========================================================
      SUBROUTINE Q2WTNEW97(Q2,WTQ2,IN)
C     ===========================================================
C
      Implicit None

#include "local.inc"
#include "common97.inc"

      LOGICAL Q2WTFIRST
      DATA  Q2WTFIRST /.TRUE./

      REAL WTQ2,WTDQ2,WTNDQ2,Q2, lumc1, lumc2, lumc3, lumc4,lumc31
      INTEGER IN

      IF (Q2WTFIRST.EQV..TRUE.) THEN
       LUMC1=LUMCND197+LUMCD197
       LUMC2=LUMCND197+LUMCD197+LUMCND297+LUMCD297
       LUMC3=LUMCND197+LUMCD197+LUMCND297+LUMCD297+LUMCND397+LUMCD397
	LUMC31=LUMCND197+LUMCD197+LUMCND297+LUMCD297+LUMCND397+LUMCD397
     &+LUMCND3197+LUMCD3197
       LUMC4=LUMCND197+LUMCD197+LUMCND297+LUMCD297+LUMCND397+LUMCD397
     & +LUMCND497+LUMCD497+LUMCND3197+LUMCD3197
       WRITE(6,*) '  *******************************'       
       WRITE(6,*) ' '       
       WRITE(6,*) ' LUMI DATA :',LUMDA97
       WRITE(6,*) ' LUMI MC :'
       WRITE(6,*) LUMC1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(6,*) LUMC2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(6,*) LUMC3,'for',Q2MC3,'< Q2 <',Q2MC31
       WRITE(6,*) LUMC31,'for',Q2MC31,'< Q2 <',Q2MC4
       WRITE(6,*) LUMC4,'for',Q2MC4,'< Q2'
       WRITE(6,*) ' '       
       WRITE(6,*) ' '       
       WRITE(6,*) ' LUMI FOR NONDIF MC :'
       WRITE(6,*) LUMCND197,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(6,*) LUMCND297,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(6,*) LUMCND397,'for',Q2MC3,'< Q2 <',Q2MC31
       WRITE(6,*) LUMCND3197,'for',Q2MC31,'< Q2 <',Q2MC4
       WRITE(6,*) LUMCND497,'for',Q2MC4,'< Q2'
       WRITE(6,*) ' '       
       WRITE(6,*) '    LUMI FOR DIF MC :'
       WRITE(6,*) LUMCD197,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(6,*) LUMCD297,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(6,*) LUMCD397,'for',Q2MC3,'< Q2 <',Q2MC31
       WRITE(6,*) LUMCD3197,'for',Q2MC31,'< Q2 <',Q2MC4
       WRITE(6,*) LUMCD497,'for',Q2MC4,'< Q2'
       WRITE(6,*) ' '       
       WRITE(6,*) ' *******************************'   
       WRITE(999,*) '  *******************************'       
       WRITE(999,*) ' '       
       WRITE(999,*) ' LUMI DATA :',LUMDA97
       WRITE(999,*) ' LUMI MC :'
       WRITE(999,*) LUMC1,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(999,*) LUMC2,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(999,*) LUMC3,'for',Q2MC3,'< Q2 <',Q2MC31
       WRITE(999,*) LUMC31,'for',Q2MC31,'< Q2 <',Q2MC4
       WRITE(999,*) LUMC4,'for',Q2MC4,'< Q2'
       WRITE(999,*) ' '       
       WRITE(999,*) ' '       
       WRITE(999,*) ' LUMI FOR NONDIF MC :'
       WRITE(999,*) LUMCND197,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(999,*) LUMCND297,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(999,*) LUMCND397,'for',Q2MC3,'< Q2 <',Q2MC31
       WRITE(999,*) LUMCND3197,'for',Q2MC31,'< Q2 <',Q2MC4
       WRITE(999,*) LUMCND497,'for',Q2MC4,'< Q2'
       WRITE(999,*) ' '       
       WRITE(999,*) '    LUMI FOR DIF MC :'
       WRITE(999,*) LUMCD197,'for',Q2MC1,'< Q2 <',Q2MC2
       WRITE(999,*) LUMCD297,'for',Q2MC2,'< Q2 <',Q2MC3
       WRITE(999,*) LUMCD397,'for',Q2MC3,'< Q2 <',Q2MC31
       WRITE(999,*) LUMCD3197,'for',Q2MC31,'< Q2 <',Q2MC4
       WRITE(999,*) LUMCD497,'for',Q2MC4,'< Q2'
       WRITE(999,*) ' '       
       WRITE(999,*) ' *******************************'   
  
	      	     
       Q2WTFIRST=.FALSE.
      ENDIF

    
      WTQ2=0.0 

c -----------------------------------
c Not Q^2 > 0.75 sample
c -----------------------------------

      if (27.52*y_tru/(27.52-hcl_eneg).le.0.04) then
      IF (Q2.GE.Q2MC1.AND.Q2.LT.Q2MC2) THEN
	 WTDQ2=LUMDA97/LUMCD197
       WTNDQ2=LUMDA97/LUMCND197

      ELSEIF (Q2.GE.Q2MC2.AND.Q2.LT.Q2MC3) THEN
       WTDQ2=LUMDA97/(LUMCD197+LUMCD297)
       WTNDQ2=LUMDA97/(LUMCND197+LUMCND297)

      ELSEIF (Q2.GE.Q2MC3.AND.Q2.LT.Q2MC4) THEN
       WTDQ2=LUMDA97/(LUMCD197+LUMCD297+LUMCD397)
       WTNDQ2=LUMDA97/(LUMCND197+LUMCND297+LUMCND397)

      ELSEIF (Q2.GE.Q2MC4) THEN
      WTDQ2=LUMDA97/(LUMCD197+LUMCD297+LUMCD397+LUMCD497)
      WTNDQ2=LUMDA97/(LUMCND197+LUMCND297+LUMCND397+LUMCND497)
      ENDIF

c ---------------------------------------------
c New bit for Q^2 > 0.75
c ---------------------------------------------

      elseif (27.52*y_tru/(27.52-hcl_eneg).gt.0.04) then
      IF (Q2.GT.Q2MC1.AND.Q2.LE.Q2MC2) THEN
	 WTDQ2=LUMDA97/LUMCD197
      WTNDQ2=LUMDA97/LUMCND197

      ELSEIF (Q2.GT.Q2MC2.AND.Q2.LE.Q2MC3) THEN
       WTDQ2=LUMDA97/(LUMCD197+LUMCD297)
       WTNDQ2=LUMDA97/(LUMCND197+LUMCND297)

      ELSEIF (Q2.GT.Q2MC2.AND.Q2.LE.Q2MC31) THEN
      WTDQ2=LUMDA97/(LUMCD197+LUMCD297+LUMCD397)
      WTNDQ2=LUMDA97/(LUMCND197+LUMCND297+LUMCND397)
 
      ELSEIF (Q2.GT.Q2MC31.AND.Q2.LE.Q2MC4) THEN
       WTDQ2=LUMDA97/(LUMCD197+LUMCD297+LUMCD397+LUMCD3197)
       WTNDQ2=LUMDA97/(LUMCND197+LUMCND297+LUMCND397+LUMCND3197)

      ELSEIF (Q2.GT.Q2MC4) THEN
      WTDQ2=LUMDA97/(LUMCD197+LUMCD297+LUMCD397+LUMCD3197+LUMCD497)
      WTNDQ2=LUMDA97/(LUMCND197+LUMCND297+LUMCND397+LUMCND3197+
     &LUMCND497)
      ENDIF
      end if
      
c ------------------------------------
c Now add diff/non diff fractions
c ------------------------------------
      IF (IN.GE.1.AND.IN.LE.NNODIFFILES97) THEN
       WTQ2=WTNDQ2*(1.0-FRACDIF97)
  
      ELSEIF (IN.GT.NNODIFFILES97.AND.IN.LE.NMCFILES97) THEN
       WTQ2=WTDQ2*FRACDIF97


      ELSEIF (IN.LT.1.OR.IN.GT.NMCFILES97) THEN
       WRITE (6,*) ' '
       WRITE (6,*) '  ERROR IN WTQ2  '
       WRITE (6,*) ' '
       WRITE (999,*) ' '
       WRITE (999,*) '  ERROR IN WTQ2  '
       WRITE (999,*) ' '

      ENDIF

      RETURN
      END                                                                        


