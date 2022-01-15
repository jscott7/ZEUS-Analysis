c ------------------------------------------------------
      subroutine count_in_bins(histofill)
c ------------------------------------------------------

      implicit none
 
#include "local.inc"
#include "constant.inc"  

	integer histofill
      real datawt, bakwt, mcwt

      datawt = 0.
	bakwt = 0.
	mcwt = 0.

	if (acceptance.ne.0) then
	   datawt = 1. / acceptance
         bakwt = (baknorm / acceptance) 
	else
	   datawt = 1.
	   bakwt = baknorm
	end if

	mcwt = mcwtq2 * f2weight * wt_vtx * bpcor

c ---------------------------
c Count data events
c ---------------------------
      if (histofill.eq.1) then
         bindat = bindat + (datawt)
	   datraw = datraw + 1

c ---------------------------
c Count background events
c ---------------------------
      else if (histofill.eq.2) then
         binbgd = binbgd + (bakwt)
	   bakraw = bakraw + 1

c ---------------------------
c Count MC events
c ---------------------------
      else if (histofill.eq.3) then
         binmc = binmc + (mcwt)
	   mcraw = mcraw + 1
	    
      end if
      
      return
      end
