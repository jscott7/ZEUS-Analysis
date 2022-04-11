c ---------------------------------------------------------
      subroutine anevent96(Filename, idnt, histofill)
c ---------------------------------------------------------
c For *** FL-ISR ANALYSIS of 1996 ntuples *** 
c ---------------------------------------------------------
      implicit none

#include"common96.inc"
#include"local.inc"

      integer histofill, calnorm
      integer istat,NEvents,i,j,ierr

      integer acceptisr, acceptempz,accept2
      integer processed, idnt

      real beamcorq,mcwt2
      real temp2

      save processed
      data processed/0/
      character*(*) Filename
      integer count
c      save count
c      data count/0/
      logical vfirst

      save vfirst
      data vfirst/.true./
      
c --- Open the data file
      write(*,*) Filename
	
      call hropen(80,'ntuple',Filename,' ',4096,istat)
      if (istat.ne.0) then
         write (*,*) 'Failed to open input file'
         stop
      endif
c
c --- Setup the common block
c
      call hrin(idnt,9999999,0)
      call hbname(idnt,' ',0,'$CLEAR')
      call hbname(idnt,' ',0,'$CLEAR')
      call hbname(idnt,flt_form(1),TrigDat,'$SET')
      call hbname(idnt,tlt_form(1),TLT,'$SET')
      call hbname(idnt,trk_form(1),VCT_XVC,'$SET')
      call hbname(idnt,cal_form(1),FEMC_EN,'$SET')
      call hbname(idnt,elec_form(1),EPROB,'$SET')
      call hbname(idnt,zufos1_form(1),ZPxBSp,'$SET')
      call hbname(idnt,zufos2_form(1),Zgamma,'$SET')
      call hbname(idnt,zufos3_form(1),ZEminPz,'$SET')
      call hbname(idnt,zufos4_form(1),ZPxCal,'$SET')
      call hbname(idnt,zufos5_form(1),ZufoE,'$SET')
      call hbname(idnt,temp_form(1),TEMPLUMG,'$SET')
      call hbname(idnt,tag1_form(1),ENE44M,'$SET')
      call hbname(idnt,tag2_form(1),TAG44_E,'$SET')
      call hbname(idnt,lumi_form(1),ENE_LE,'$SET')
      call hbname(idnt,lumi2_form(1),ENE_LG,'$SET')
      call hbname(idnt,gen_form(1),RUN_NUM,'$SET')
      call hbname(idnt,mctrue_form(1),MCELE,'$SET')
      call hbname(idnt,bgdtup_form(1),EVTwant,'$SET')
      
      call hnoent(idnt,NEvents)

      write (*,*) 'ANALYSE: File contains ',NEvents,' Events'
      write (999,*) 'ANALYSE: File contains ',NEvents,' Events'

      count = 0
c ----------------------------------------
c Loop over contents of this ntuple
c ----------------------------------------    
      do i=1,NEvents

c ----------------------------------------
c Fill common block with event data
c ----------------------------------------
         if (i.eq.1) then
            call hgnt(idnt,i,ierr)
         else
            call hgntf(idnt,i,ierr)
         endif

         if (ierr.ne.0) then
            write (*,*) 'ANALYSE: Failed to get event ',i
            write (999,*) 'ANALYSE: Failed to get event ',i
         else

c ---------------------------------------------
c 3 bad events in MC (undefined CAL energy) - Reject
c ---------------------------------------------
	 if ((cal_xp**2+cal_yp**2).lt.0.) then
            write(*,*) 'ERROR: NTUPLE INPUT BAD'
	    goto 10
	 end if

c ------------------------------------------
c Do ISR background addition
c ------------------------------------------

         call bgdlumi96(histofill)

         WT_VTX = 1.0

	 if (histofill.eq.3) then

C ===========================================
C === GET A WEIGHT FOR VERTEX REWEIGHTING ===
C ===========================================

c               IF (YEAR.EQ.1996) THEN
c                CALL VTX_WEIGHT(MC_ZV,WT_VTX,963,96,0)
            call Rewt_Vtx_Isr(MC_ZV,0,WT_VTX,Ierr)
c                ELSEIF (YEAR.EQ.1997) THEN
c                 CALL VTX_WEIGHT(MC_ZV,WT_VTX,971,97,0)
c               ENDIF

c ------------------------------------------
c Do lumi correction for ISR MC
c ------------------------------------------            

            call corrmclumig96         
         end if

c ------------------------------------------
c General corrections and calculations F2ISR
c ------------------------------------------
         call corrections96
         electron_en = corrected_en
         electron_th = best_th

         call isr_calculate96

c ------------------------------------------
c Now apply cuts
c ------------------------------------------               
         acceptisr = 0
	 acceptempz = 0
	 calnorm = 0
         accept2 = 0

         call isr_cuts96(acceptisr, acceptempz, accept2, calnorm,
     &                    histofill)

c ------------------------------------------
c Calculate Weights
c ------------------------------------------         
         if ((acceptisr.eq.1).or.(acceptempz.eq.1)) then

c ------------------------------------------
c Montecarlo
c ------------------------------------------
            if (histofill.eq.3) then
	       if (year.eq.1996) then
                  call Q2WTNEW96(q2_tru, mcwtq2, mcinput)
	       end if

               call f2wt96
	       bpcor = beamcorq(Zgamma,q2_el,2,ierr)
c	        bpcor = 1.0
            else

c ------------------------------------------
c Data / Background
c ------------------------------------------
               call lumiacceptance(RUN_NUM, acceptance)
            end if
	 end if

c ----------------------------------------------
c Normalise background - for total E-Pz > 62GeV.
c ----------------------------------------------
	 if (calnorm.eq.1) then
            if (empztot.gt.bgdnormal) then
	       if (histofill.eq.1) then
                  datanorm = datanorm + 1.
	       else if (histofill.eq.20) then
		  backnorm = backnorm + 1.
	       end if
	    end if
	 end if

c ------------------------------------------
c Fill histos & count in bins
c ------------------------------------------

c ------------------------------------------
c For all events 
c ------------------------------------------
         if (accept2.eq.1) then
	    if (histofill.eq.3) then
	       temp2 = (ysigz-y_tru)/y_tru
	       mcwt2 = mcwtq2 * f2weight * wt_vtx * bpcor
	      
               if (log10(y_tru).gt.-1.7
     &              .and.log10(y_tru).le.-1.55) then
                 call hf1(430,temp2,mcwt2)
               elseif (log10(y_tru).gt.-1.55
     &              .and.log10(y_tru).le.-1.4) then
                 call hf1(431,temp2,mcwt2)
               elseif (log10(y_tru).gt.-1.4
     &              .and.log10(y_tru).le.-1.25) then
                 call hf1(432,temp2,mcwt2)
               elseif (log10(y_tru).gt.-1.25
     &              .and.log10(y_tru).le.-1.1) then
	         call hf1(433,temp2,mcwt2)
               elseif (log10(y_tru).gt.-1.1
     &              .and.log10(y_tru).le.-0.95) then
	         call hf1(434,temp2,mcwt2)
               elseif (log10(y_tru).gt.-0.95
     &              .and.log10(y_tru).le.-0.8) then
	         call hf1(435,temp2,mcwt2)
               elseif (log10(y_tru).gt.-0.8
     &              .and.log10(y_tru).le.-0.65) then
	         call hf1(436,temp2,mcwt2)
               elseif (log10(y_tru).gt.-0.65
     &              .and.log10(y_tru).le.-0.5) then
	         call hf1(437,temp2,mcwt2)
               elseif (log10(y_tru).gt.-0.5
     &              .and.log10(y_tru).le.-0.35) then
	         call hf1(438,temp2,mcwt2)
               elseif (log10(y_tru).gt.-0.35
     &              .and.log10(y_tru).le.-0.2) then
	         call hf1(439,temp2,mcwt2)
               end if
	    end if
         end if

c ------------------------------------------
c For events within FL bin
c ------------------------------------------
         if (acceptisr.eq.1) then    
	    call count_in_bins(histofill)
            call isr_histos96(histofill)
            count = count + 1    
	 end if
           
c ------------------------------------------
c Fill E-Pz Histos
c ------------------------------------------
         if (acceptempz.eq.1) then
            call empzfill96(histofill)
         end if

      end if

10	continue
      enddo
    
      write(*,*) 'ISR events', count
      write(999,*) 'ISR events', count

c ------------------------------------------     
c All done, so clean up this file
c ------------------------------------------    
      call hrend('ntuple')

      close(80)
     
      return
      end
