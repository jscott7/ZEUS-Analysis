c ---------------------------------------------------
      subroutine isr_cuts96(acceptisr, acceptempz, accept2,
     &                    calnorm, histofill)
c ---------------------------------------------------
c *** FL CUTS ***
c ---------------------------------------------------
      implicit none
#include "local.inc"
#include "common96.inc"
#include "constant.inc"


      integer acceptisr, acceptempz, calnorm,accept2
      integer histofill
      logical boxcut,h1box
      integer elumiecut, elumigcut
      logical offlinecut
      integer probcut, empzcalcut, vtxcut, empztotcut
      integer energycut
      integer yjbcut, ycut, yjbcorrcut, ycorrcut
      integer sltempzcut
      integer tagger44mcut
      integer q2cut
      integer deltaisrcut
      integer etamaxcut
      integer fltflag
      integer bit, tltbit
      logical tltcut
      integer bincut
      integer lumigeomcut
      integer haccut,bgdcut,hacencut
      real boxadd,zufoesc
      real temp

      logical firstisrc
      data firstisrc/.true./
      save firstisrc

      acceptisr = 0
      acceptempz = 0
      calnorm = 0
      haccut = 0
      hacencut = 0
      bgdcut = 0
      lumigeomcut = 0
      boxcut = .false.
      elumiecut = 0
      elumigcut = 0
      tltcut = .false.
      probcut = 0
      empzcalcut = 0
      empztotcut = 0
      vtxcut = 0
      energycut = 0
      yjbcut = 0
      ycut = 0
      ycorrcut = 0
      yjbcorrcut = 0
      sltempzcut = 0
      tagger44mcut = 0
      q2cut = 0
      deltaisrcut = 0
      etamaxcut = 0
      bincut = 0

      if (firstisrc) then
         write(*,*) '***************************'
         write(*,*) '***  Start FL/ISR cuts  ***'
         write(*,*) '***************************'
         write(999,*) '***************************'
         write(999,*) '***  Start FL/ISR cuts  ***'
         write(999,*) '***************************'
         firstisrc = .false.
      end if     

c --------------------------------
c 13 x 7 Box Cut if SRTD available
c else 15 x 15 Box Cut
c 
c From Adi
c -------------------------------

      boxcut = h1box(bestx,besty)

c -------------------------------
c Outer boxcut 30 x 30 cm 
c Essentially only use SRTD hit.
c -------------------------------

      if ((abs(bestx).gt.30.).or.(abs(besty).gt.30.)) then
         boxcut = .false.
      end if   

c -------------------------------
c Lumi Cuts
c -------------------------------
      if ((elumie.lt.2.).and.(elumie.gt.-1000.)) then
         elumiecut = 1
      end if

      if ((elumig.gt.elglow).and.(elumig.lt.30.)) then
         elumigcut = 1
      end if

c ----------------------------------------
c Use offline trigger only for Data and MC
c Background forced to true. (DIS02)
c ----------------------------------------
	
      if (offlineuse.eq.1) then
      if ((histofill.eq.2).or.(histofill.eq.20)) then
      if (iand(tlt(4),131072).eq.131072) then	
         tltcut = .true.
      else	
         tltcut =.false.
      end if

      else
         call triggercut96(tltcut, fltflag)
      end if

      else
c --------------------------------------
c For background require DIS02, else DIS10
c --------------------------------------

         bit = 10
	 tltbit = 2**(bit-1+16)

c Background DIS02 
	if ((histofill.eq.2).or.(histofill.eq.20)) then
         if (iand(tlt(4),131072).eq.131072) then
	    tltcut = .true.
	 else
	    tltcut = .false.
	 end if
c Data DIS10
      else if (histofill.eq.1) then
          if (iand(tlt(4),tltbit).eq.tltbit) then
	      tltcut = .true.
	  end if
c MC DIS10
c Dont apply TLT cut for MC as SRTD_GOOD info not simulated.
      else if (histofill.eq.3) then
c	    if (iand(tlt(4),tltbit).eq.tltbit) then
		  tltcut = .true.
c           end if
      else
         tltcut = .false.
      end if

      end if

c ------------------------------------
c E - Pz Cuts
c E - Pz Cal > 22 GeV for Fl analysis
c ------------------------------------
      if (empzcal.gt.empzcalcard) then
         empzcalcut = 1
      end if

c ---------------------------------
c Total E-Pz => 48 - 60 GeV for '96
c               46 - 60 GeV for '97
c Now from steering cards.
c ---------------------------------

      if ((empztot.gt.empztotlow).and.
     &    (empztot.lt.empztothigh)) then
        empztotcut = 1
      end if

      if (trigdat(4).gt.0.0) then
         sltempzcut = 1
      end if
	
c -------------------------------
c y cuts (use Zufos for y_jb)
c Changes for Fl analysis : 
c y_jb > 0.0
c y_jbcorr > 0.05
c y_elcorr > 0.05 
c -------------------------------
      if ((y_el.gt.-999.).and.(y_el.lt.0.95)) then
         ycut = 1
      end if

      if ((y_elcorr.gt.0.05).and.(y_elcorr.lt.0.95)) then
         ycorrcut = 1
      end if

      if ((zy_jb.gt.0.0).and.(zy_jb.lt.0.95)) then
         yjbcut = 1
      end if

      if ((zy_jbcorr.gt.0.05).and.(zy_jbcorr.lt.0.95)) then
         yjbcorrcut = 1
      end if

c -------------------------------
c Others 
c -------------------------------

      if ((EPROB.gt.0.9).and.(EPROB.le.1.0)) then
         probcut = 1
      end if

      if (abs(VCT_ZVC).lt.zvertex) then
         vtxcut = 1
      end if

      if ((electron_en.gt.enelow).and.(electron_en.lt.30.)) then
         energycut = 1
      end if

      if (tag44mveto.eq.1) then
         if (en44m.lt.60.) then
            tagger44mcut = 1
         end if
      else 
         tagger44mcut = 1
      end if

      if ((q2corr.gt.0.1).and.(q2corr.lt.90000.0)) then
         q2cut = 1
      end if

      if ((deltaisr.gt.-10.).and.(deltaisr.lt.10.)) then
         deltaisrcut = 1
      end if

      if ((ETAMAX.gt.-10.).and.(ETAMAX.lt.1000.)) then
         etamaxcut = 1
      end if

      if ((q2corr.ge.q2min).and.(q2corr.lt.q2max).and.(ysigz.ge.ymin)
     &     .and.(ysigz.lt.ymax)) then
         bincut = 1
      end if

c ---------------------------------
c Lumi geometric acceptance cut.
c For True MC
c ---------------------------------
      if (histofill.eq.3) then
	 call lumigeom96(lumigeomcut)
      else
         lumigeomcut = 1
      end if

      if (HACenrgy.lt.5) haccut = 1

      zufoesc=zufoe*hadscale
      if (zufoesc.gt.2) hacencut = 1

      if ((qedc_no.eq.0).and.(muon_no.eq.0)) then
         bgdcut = 1
      end if

c -------------------------------
c Apply all the cuts
c -------------------------------

      if ( 
     &	 (tltcut).and.
     &     (lumigeomcut.eq.1).and.
     &     (boxcut).and.
     &     (Elumiecut.eq.1).and.
     &     (Elumigcut.eq.1).and.
     &     (probcut.eq.1).and.
     &     (vtxcut.eq.1).and.
     &     (empzcalcut.eq.1).and.
c ---------------------------------
     &     (empztotcut.eq.1).and.
c ---------------------------------
     &     (sltempzcut.eq.1).and.
c     &     (yjbcut.eq.1).and.
     &     (yjbcorrcut.eq.1).and.
c     &     (ycut.eq.1).and.
     &	 (ycorrcut.eq.1).and.
     &     (energycut.eq.1).and.
c     &     (deltaisrcut.eq.1).and.
c     &	 (etamaxcut.eq.1).and.
     &	 (q2cut.eq.1).and.
     &     (tagger44mcut.eq.1).and.
     &     (haccut.eq.1).and.
     &     (hacencut.eq.1).and.
     &     (bgdcut.eq.1).and.
     &     (bincut.eq.1)
     &     ) then

         acceptisr = 1

      end if

c --------------------------------
c Fill  E-Pz plots
c --------------------------------

      if ( 
     &     (tltcut).and.
     &     (lumigeomcut.eq.1).and.
     &     (boxcut).and.
     &     (Elumiecut.eq.1).and.
     &     (Elumigcut.eq.1).and.
     &     (probcut.eq.1).and.
     &     (vtxcut.eq.1).and.
     &     (empzcalcut.eq.1).and.
     &     (sltempzcut.eq.1).and.
c     &     (yjbcut.eq.1).and.
     &     (yjbcorrcut.eq.1).and.
c     &     (ycut.eq.1).and.
     &	 (ycorrcut.eq.1).and.
     &     (energycut.eq.1).and.
c     &     (deltaisrcut.eq.1).and.
c     &	 (etamaxcut.eq.1) .and.
     &	 (q2cut.eq.1).and.
     &     (tagger44mcut.eq.1).and.
     &     (haccut.eq.1).and.
     &     (hacencut.eq.1).and.
     &     (bgdcut.eq.1).and.
     &     (bincut.eq.1) 
     &     ) then
      
	   acceptempz = 1

	   if ((histofill.eq.20).or.(histofill.eq.1)) then
	      calnorm = 1
	   end if

	end if

c ------------------------------------------
c Whole range plots
c ------------------------------------------
      if ( 
     &	 (tltcut).and.
     &     (lumigeomcut.eq.1).and.
     &     (boxcut).and.
     &     (Elumiecut.eq.1).and.
     &     (Elumigcut.eq.1).and.
     &     (probcut.eq.1).and.
     &     (vtxcut.eq.1).and.
     &     (empzcalcut.eq.1).and.
c ---------------------------------
     &     (empztotcut.eq.1).and.
c ---------------------------------
     &     (sltempzcut.eq.1).and.
c     &     (yjbcut.eq.1).and.
     &     (yjbcorrcut.eq.1).and.
c     &     (ycut.eq.1).and.
     &	 (ycorrcut.eq.1).and.
     &     (energycut.eq.1).and.
c     &     (deltaisrcut.eq.1).and.
c     &	 (etamaxcut.eq.1).and.
     &	 (q2cut.eq.1).and.
     &     (tagger44mcut.eq.1).and.
     &     (haccut.eq.1).and.
     &     (hacencut.eq.1).and.
     &     (bgdcut.eq.1)
     &     ) then

         accept2 = 1

      end if

c -----------------------------------------
c Kinematic peak
c ------------------------------------------

      if (
     & (boxcut).and.
     & (probcut.eq.1).and.
     & (zy_jbcorr.lt.0.03).and.
     & (zy_jbcorr.gt.0.01)
     & )then

	temp = electron_en + elumig
	if (histofill.eq.1) then
	if ((electron_en.ge.4.).and.(electron_en.lt.7)) then
	call hf1(500,temp,1.)
	else if((electron_en.ge.7.).and.(electron_en.lt.10)) then
	call hf1(501,temp,1.)
	else if((electron_en.ge.10.).and.(electron_en.lt.13)) then
	call hf1(502,temp,1.)
	else if((electron_en.ge.13.).and.(electron_en.lt.16)) then
	call hf1(503,temp,1.)
	else if((electron_en.ge.16.).and.(electron_en.lt.19)) then
	call hf1(504,temp,1.)
	else if((electron_en.ge.19.).and.(electron_en.lt.22)) then
	call hf1(505,temp,1.)
	end if

	else if (histofill.eq.3) then
	if ((electron_en.ge.4.).and.(electron_en.lt.7)) then
	call hf1(506,temp,1.)
	else if((electron_en.ge.7.).and.(electron_en.lt.10)) then
	call hf1(507,temp,1.)
	else if((electron_en.ge.10.).and.(electron_en.lt.13)) then
	call hf1(508,temp,1.)
	else if((electron_en.ge.13.).and.(electron_en.lt.16)) then
	call hf1(509,temp,1.)
	else if((electron_en.ge.16.).and.(electron_en.lt.19)) then
	call hf1(510,temp,1.)
	else if((electron_en.ge.19.).and.(electron_en.lt.22)) then
	call hf1(511,temp,1.)
	end if

	end if

      end if

      return
      end
	

