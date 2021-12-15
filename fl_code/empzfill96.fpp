c ---------------------------------------
      Subroutine empzfill96(histofill)
c ---------------------------------------

      implicit none

#include "local.inc"
#include "common96.inc"
#include "constant.inc"

      integer histofill
	real datawt, bakwt, mcwt

	if (acceptance.ne.0.) then
         datawt = 1. / acceptance
         bakwt = baknorm / acceptance
	else
	   datawt = 0.
	   bakwt = 0.
	end if

	mcwt = mcwtq2 * f2weight * wt_vtx * bpcor

c ---------------------------------        
c  Data
c ---------------------------------
         if (histofill.eq.1) then
            call hf1(103,zempzhad,datawt)
            call hf1(104,empzcal,datawt)
            call hf1(105,empztot,datawt)
            
c ------------------------------------
c Background
c ------------------------------------
         else if (histofill.eq.2) then
            call hf1(203,zempzhad,bakwt)
            call hf1(204,empzcal,bakwt)
            call hf1(205,empztot,bakwt) 
 
c ------------------------------------   
c MonteCarlo
c ------------------------------------        
         else if (histofill.eq.3) then
            call hf1(303,zempzhad,mcwt)
            call hf1(304,empzcal,mcwt)
            call hf1(305,empztot,mcwt)            
         end if

      return
      end




