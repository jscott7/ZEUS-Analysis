C
C    ============================ 
      Subroutine LUMISIM(SIMMODE)
C    ============================ 

      Implicit None

#include "common.inc"
#include "local.inc"

      real off,nl,res,off1,offold,nlold,resold,
     + off96,nl96,res96,off97,nl97,res97,
     + off97n,nl97n,res97n,resad97n,alpha97n,
     + alpha,beta,resad
	integer SIMMODE
 
      data off96,nl96,res96 /0.2,0.0011,0.23/
 
c      data off97n,nl97n,res97n,resad97n,alpha97n 
c     + /0.38,-0.0005,0.323,0.3,19.5/

      real Etrue,arnd,brnd,xgam
      real rngama
	real rangam
	real Eres
      real Em1,em2,em3,off2
	real DEF_EE
	data DEF_EE /27.55/

      Logical LumiSimFirst
      Logical LSFMC
      Logical LSFBGD

      Data LumiSimFirst /.true./
      Data LSFMC /.true./
      Data LSFBGD /.true./

      SAVE LumiSimFirst,LSFMC,LSFBGD
	off97n = 0.38
	nl97n = -0.0005
	res97n = 0.323
	resad97n = 0.3
	alpha97n = 19.5

      IF (SIMMODE.EQ.1) THEN
       NL=NL96
       OFF=OFF96
       RES=RES96
 
      ELSEIF (SIMMODE.EQ.3) THEN
       NL=NL97N
       OFF=OFF97N
       RES=RES97N
       RESAD=RESAD97N
       ALPHA=ALPHA97N
      ENDIF

      If(LumiSimfirst)then
         write(6,*) '  '
         WRITE(6,*) ' ***************************************'
         write(6,*) '  '
         write(6,*) '      First call of LUMI SIMULATION '
         write(6,*) '           and LUMISMEARING'
         write(6,*) '  '
         write(6,*) '         RESOLUTION :', RES
         write(6,*) '             OFFSET :', OFF
         write(6,*) '       NONLINEARITY :', NL
         write(6,*) ' EXTRA RESOLUTION 1 :', RESAD
         write(6,*) '          ASYMMETRY :', ALPHA 
         write(6,*) '  '
         write(6,*) '        LUMISIMTYP :',SIMMODE
         write(6,*) '       BEAM ENERGY :',DEF_EE
         write(6,*) '  '
         write(999,*) '  '
         WRITE(999,*) ' ***************************************'
         write(999,*) '  '
         write(999,*) '      First call of LUMI SIMULATION '
         write(999,*) '           and LUMISMEARING'
         write(999,*) '  '
         write(999,*) '         RESOLUTION :', RES
         write(999,*) '             OFFSET :', OFF
         write(999,*) '       NONLINEARITY :', NL
         write(999,*) ' EXTRA RESOLUTION 1 :', RESAD
         write(999,*) '          ASYMMETRY :', ALPHA 
         write(999,*) '  '
         write(999,*) '        LUMISIMTYP :',SIMMODE
         write(999,*) '       BEAM ENERGY :',DEF_EE
         write(999,*) '  '         
         
	   LumiSimfirst =.false.
      endif

C call random generator
C
      call rannor(arnd,brnd)
c	call rnorml(brnd,1)
C
C
      IF (SIMMODE.EQ.1) THEN
       If(LSFMC)then
         WRITE(6,*) ' ***************************************'
         WRITE(6,*) ' '        
         WRITE(6,*) '    FIRST CALL TO LUMISIM IN MC MODE'        
         WRITE(6,*) ' '        
         WRITE(6,*) ' ***************************************'
         LSFMC =.false.
       EndIf
	end if

	elumig = 0.
	Etrue = 0.
       
	Etrue = HCL_ENEG
      if(etrue.le.0.)  then
         elumig = 0.
         return
      endif

c ---------------------------------------
c 1996
c ---------------------------------------
       IF (SIMMODE.EQ.1) THEN
        em1=(Etrue-off)*(1.0+nl*(DEF_EE-Etrue))
        elumig=em1 + res*sqrt(Etrue)*brnd
       ENDIF

c ---------------------------------------
c 1997
c ---------------------------------------	  
c	beta = 0.
c	xgam = 0.
c	em1 = 0.
c	em2 = 0.
       IF (SIMMODE.EQ.3) THEN

      etrue = etrue - off
 	if (etrue.gt.0.0) then
	 beta = sqrt(alpha/(etrue*res**2+resad))
	 xgam = rngama(alpha)/beta
	 etrue = etrue - xgam+(alpha-1.)/beta
	 if (etrue.lt.0) etrue = 0.
	end if						
      elumig = etrue*(1.0+nl*(DEF_EE-off-etrue))

c       beta=sqrt(alpha/(((ETRUE-OFF)*res**2)+resad))
c       xgam=(rngama(alpha))/beta
c	 em1 = (ETRUE-OFF-xgam+(alpha-1.0)/beta)
c	 em2 = (1.0+nl*(DEF_EE-ETRUE+xgam-(alpha-1.0)/beta))
c     elumig= em1 * em2
	
c       em1=(Etrue-off)*(1.0+nl*(DEF_EE-Etrue))
c       elumig=em1 + res*sqrt(Etrue)*brnd
       ENDIF

      return
      end


