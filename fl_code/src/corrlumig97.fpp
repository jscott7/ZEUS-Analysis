c ----------------------------
      subroutine corrmclumig97
c ----------------------------
c Simulate the energy response of Lumig - simulation in MC is wrong.

      implicit none
      logical firstc
      data firstc/.true./
      save firstc

#include"common97.inc"
#include"local.inc"
#include"constant.inc"

      real off,nl,res,off1,offold,nlold,resold,
     + off96,nl96,res96,off97,nl97,res97,
     + off97n,nl97n,res97n,resad97n,alpha97n,
     + alpha,beta,resad

      data off96,nl96,res96 /0.45,0.0011,0.32/
 
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
      real x

      if (firstc) then
         write(*,*) '********************************'
         write(*,*) '***  Lumi correction for MC  ***'
         write(*,*) '********************************'
         write(999,*) '********************************'
         write(999,*) '***  Lumi correction for MC  ***'
         write(999,*) '********************************'
	   write(*,*) 'YEAR SET TO:',year
         firstc = .false.
	end if

	elumie = 0.
	elumig = 0.
	en44m = 0.

	off97n = 0.38
	nl97n = -0.0005
	res97n = 0.323 + lumires
	resad97n = 0.3
	alpha97n = 19.5

c      IF (year.eq.1996) THEN
c       NL=NL96
c       OFF=OFF96
c       RES=RES96
 
c      ELSEIF (year.eq.1997) THEN
       NL=NL97N
       OFF=OFF97N
       RES=RES97N
       RESAD=RESAD97N
       ALPHA=ALPHA97N
c     ENDIF

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
         write(6,*) '        LUMISIMTYP :',year
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
         write(999,*) '        LUMISIMTYP :',year
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
c       IF (year.EQ.1996) THEN
c        em1=(Etrue-off)*(1.0+nl*(DEF_EE-Etrue))
c        elumig=em1 + (res)*sqrt(Etrue)*brnd
c       ENDIF
c	elumig = elumig + lumiscale

c ---------------------------------------
c 1997
c ---------------------------------------	  
	beta = 0.
	xgam = 0.
	em1 = 0.
	em2 = 0.
c     IF (year.EQ.1997) THEN
c
c       etrue = etrue - off
c	 if (etrue.gt.0.0) then
c	  beta = sqrt(alpha/(etrue*res**2+resad))
c	  xgam = rngama(alpha)/beta
c	  etrue = etrue - xgam+(alpha-1.)/beta
c	  if (etrue.lt.0) etrue = 0.
c	 end if						
c     elumig = etrue*(1.0+nl*(DEF_EE-off-etrue))

      beta=sqrt(alpha/(((ETRUE-OFF)*res**2)+resad))
      xgam=(rngama(alpha))/beta
	em1 = (ETRUE-OFF-xgam+((alpha-1.0)/beta))
	em2 = (1.0+nl*(DEF_EE-ETRUE+xgam-((alpha-1.0)/beta)))
      elumig= (em1 * em2) + lumiscale
	
c      ENDIF

      elumie = ene_le
      en44m = TAG44_E      
 
c -------------------------------
c Log MC quantities
c -------------------------------
      if (Q2_TRU.gt.0.0) then
         logmcq2 = log10(Q2_TRU)
      else
         logmcq2 = -99.
      end if
         
      if (X_TRU.gt.0.0) then 
         logmcx = log10(X_TRU)
      else
         logmcx = -99.
      end if

      if (Y_TRU.gt.0.0) then
         logmcy = log10(Y_TRU)
      else
         logmcy = -99.
      end if

      return
      end 
