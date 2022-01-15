c ----------------------------------
      subroutine isr_calculate96
c ----------------------------------

      implicit none

#include "local.inc"
#include "common96.inc"
#include "constant.inc"

      real totcalen, bottom
      real temp
      real s
	real xcorr
	real def_ee, def_ep, eg_cal
	real electron_pz
	real pt2
	integer i, j
	real epzlocal

	def_ee = 0.
	eg_cal = 0.
	def_ee = 27.52
	def_ep = 820.0
	z = 0.
	Q2_el = 0.
	y_el = 0.
	y_elcorr = 0.
	zy_jb = 0.
	zy_jbcorr = 0.
	y_sig = 0.
	electron_pz = 0.
	q2corr = 0.
	deltaisr = 0.
	empzcal = 0.
	electron_pz = 0.

c ---------------------------------
c Zufos e-pz hadronic
c ---------------------------------
	if (bkspuse.eq.1) then
         zempzhad = ZEminPz
	else
	   zempzhad = (ZufoE*hadscale)-ZufoPz
	end if

c ---------------------------------
c Calculate Cal e-pz
c ---------------------------------
	electron_pz = electron_en * cos(electron_th)
	empzcal = zempzhad + electron_en - electron_pz  

c ---------------------------------
c Calculate total e-pz
c ---------------------------------
      empztot = empzcal + 2.*(elumig)

c ---------------------------------
c Calculate Q2 with electron method
c ---------------------------------
      Q2_el = 2. * def_ee * electron_en * 
     &      (1. + cos(electron_th))

c ---------------------------------
c Calculate y jaquet-blondel (Zufos)
c ---------------------------------
      zy_jb = zempzhad / (2. * def_ee)

c ---------------------------------
c Calculate y with electron method
c ---------------------------------
      y_el = 1.-((electron_en/(2.*def_ee))*
     &	(1.-cos(electron_th)))
    
c ---------------------------------
c Calculate y with sigma method
c ---------------------------------
      if ((zempzhad + (electron_en * (1. - cos(electron_th)))).ne.0.)
     & then
      y_sig = zempzhad / (zempzhad + 
     &          (electron_en * (1. - cos(electron_th))))
      end if

c -----------------------------------------
c Correct y_el for ISR (z correction factor)
c -----------------------------------------
      z = (def_ee-(elumig)) / def_ee
 
	ysigz = y_sig * z

c ---------------------------------
c y_el corr = (y_el + z -1) / z
c ---------------------------------
      if (z.ne.0.) then
      y_elcorr = (y_el + z - 1.) / z
      end if      

c ---------------------------------
c y_jb corr = y_jb / z
c ---------------------------------
      if (z.ne.0.) then
      zy_jbcorr = zy_jb / z
      end if
  
c ---------------------------------
c Correct Q^2 for ISR = zQ^2
c ---------------------------------
      q2corr = z * Q2_el

c ---------------------------------
c Calculate x = Q2 / sy
c ---------------------------------
      s = 4. * def_ee * ep
      
      temp = s * y_el
      if (temp.ne.0.) then
         x_el = Q2_el / (temp)
      end if

c --------------------------------
c Correct x for ISR
c --------------------------------
	if ((y_el+z-1.).gt.0.0) then
	   xcorr = (z * x_el * y_el) / (z + y_el - 1)
	end if

c --------------------------------
c Calculate delta_isr
c --------------------------------
	eg_cal = def_ee*(y_el - y_sig)
	if (elumig.ne.0.) then
	deltaisr = (elumig-eg_cal)/(elumig)
      end if

c --------------------------------
c Calculate Gamma Had, no backsplash
c I believe Z gamma includes backsplash....
c --------------------------------
      gamma = Zgamma

	pt2 = (zufopx**2)+(zufopy**2)
	epzlocal = zufoE-ZufoPz
	if ((pt2+epzlocal**2).ne.0) then
	gamma2=ACOS((pt2-epzlocal**2)/(pt2+epzlocal**2))
	end if
	gamma = gamma * 180./pi
	gamma2 = gamma2 * 180./pi

c --------------------------------
c Logarithmic quantities
c --------------------------------
	if (ysigz.gt.0.0) then
	   ysig_x_z = log10(ysigz)
	else
	   ysig_x_z = -99.
	end if

      if (y_elcorr.gt.0.0) then
         logyelcorr = log10(y_elcorr)
      else
         logyelcorr = -99.
      end if
      
      if (y_el.gt.0.0) then
         logyel = log10(y_el)
      else
         logyel = -99.
      end if
      
      if (y_sig.gt.0.0) then
         logysig = log10(y_sig)
      else
         logysig = -99.
      end if
      
      if (zy_jbcorr.gt.0.0) then
         logyjbz = log10(zy_jbcorr)
      else
         logyjbz = -99.
      end if

      if (xcorr.gt.0.0) then
         logxel = log10(xcorr)
      else
         logxel = -99.
      end if

      if (q2corr.gt.0.0) then
         logq2elc = log10(q2corr)
      else
         logq2elc = -99.
      end if

	if (q2_el.gt.0.0) then
         logq2el = log10(q2_el)
      else
         logq2el = -99.
      end if

      return
      end
