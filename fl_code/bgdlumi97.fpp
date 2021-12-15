c -----------------------------------------
      subroutine bgdlumi97(histofill)
c -----------------------------------------

      implicit none

#include"common97.inc"
#include"local.inc"

      integer histofill
      logical firstb 
      data firstb/.true./
      save firstb

      if (firstb) then
         if ((histofill.eq.2).or.(histofill.eq.20)) then
            write(*,*) ' **************************************'
            write(*,*) ' ***   Brems and BGD lumi addition  ***'
            write(*,*) ' **************************************'
            write(999,*) ' **************************************'
            write(999,*) ' ***   Brems and BGD lumi addition  ***'
            write(999,*) ' **************************************'
         end if
         firstb = .false.
      end if

c -----------------------------------
c Add the Brems and bgd lumi into
c new variables, locallumie, locallumig
c -----------------------------------
c 1. If data use values from ntuple.
c Multiply Lgamma by 0.36% for '96 data
c and by 1% for '97 data
c -----------------------------------
	if (histofill.eq.1) then
c	   if (year.eq.1996) then
c	      elumig = ENE_LG * 1.0036
c	   else if (year.eq.1997) then
		  elumig = ENE_LG * 1.01
c	   end if
         elumie = ENE_LE
	   en44m = TAG44_E
	end if

c -------------------------------------
c If Background add the brems lumi.
c Do same Lgamma factor.
c -------------------------------------
      if ((histofill.eq.2).or.(histofill.eq.20)) then
c	   if (year.eq.1996) then
c	      elumig = (ENE_LG + TEMPLUMG) * 1.0036
c	   else if (year.eq.1997) then
	      elumig = (ENE_LG + TEMPLUMG) * 1.01
c	   end if
         elumie = ENE_LE + TEMPLUME
	   en44m = TAG44_E + ene44m	
      end if

      return
      end
