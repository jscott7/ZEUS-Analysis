c -------------------------------------
      subroutine calculate_f2
c -------------------------------------

      implicit none

#include "local.inc"

      integer i, j
      logical first
      real a, b, c
      real f2_allm
      real err_data, err_bak, err_mc

      f2isr = 0.
      f2_allm = 0.
      a = 0.
      b = 0.
      c = 0.
      errors = 0.
      err_data = 0.
      err_bak = 0.
      err_mc = 0.
	 
      call f2allm(xmean, q2mean, f2_allm)

c -----------------------------------------
c Calculate F2
c -----------------------------------------
      f2isr = (bindat / (binmc + binbgd)) * f2_allm
      write(*,*) 'F2:', f2isr
      write(*,*) 'dat:',bindat
      write(*,*) 'bak:',binbgd
      write(*,*) 'mc:',binmc
      write(*,*) 'f2_allm:',f2_allm
      write(999,*) 'F2:', f2isr
      write(999,*) 'dat:',bindat
      write(999,*) 'bak:',binbgd
      write(999,*) 'mc:',binmc
      write(999,*) 'f2_allm:',f2_allm

c -----------------------------------------
c Calculate Statistical errors 
c -----------------------------------------
      if (datraw.gt.0.) then
         err_data = (sqrt(datraw)) * (bindat/datraw)
      else
         err_data = 0.
      end if

      if (bakraw.gt.0.) then
         err_bak = (sqrt(bakraw)) * (binbgd/bakraw)
      else
         err_bak = 0.
      end if

      if (mcraw.gt.0.) then
         err_mc = (sqrt(mcraw)) * (binmc/mcraw)
      else
         err_mc = 0.
      end if

      if ((bindat+binbgd).gt.0) then
         a = ((err_data + err_bak) / (bindat + binbgd))**2
      else 
         a = 0.
      end if

      if (binmc.gt.0) then
         b = (err_mc / binmc)**2
      else
         b = 0.
      end if

      c = sqrt(a + b)
      errors = c * f2isr
       
      write(998,*) f2isr,errors  
      return
      end

