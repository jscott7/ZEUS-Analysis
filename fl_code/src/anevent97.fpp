c ---------------------------------------------------------
      subroutine anevent97(Filename, idnt, histofill)
c ---------------------------------------------------------
c For *** FLISR ANALYSIS *** 
c ---------------------------------------------------------
      implicit none

#include"common97.inc"
#include"local.inc"

      integer histofill, calnorm
      integer istat,NEvents,i,j,ierr
	integer acceptisr, acceptempz,kpaccept
      integer processed, idnt
	real beamcorq
	real temp,datawt,mcwt

	real l(4), lp(4), p(4), v1(4), ans

      save processed
      data processed/0/
      character*(*) Filename
      integer count
c      save count
c      data count/0/
	logical vfirst
	save vfirst
	data vfirst/.true./
	real xmin,xmax,ymean

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

      call hbname(idnt,flt_form(1),TrigDat,'$SET')
      call hbname(idnt,tlt_form(1),TLT,'$SET')
      call hbname(idnt,trk_form(1),VCT_XVC,'$SET')
      call hbname(idnt,cal_form(1),FEMC_EN,'$SET')
      call hbname(idnt,elec_form(1),EPROB,'$SET')
      call hbname(idnt,zufos_form(1),ZPxBSp,'$SET')
      call hbname(idnt,lumi_form(1),TEMPLUMG,'$SET')
      call hbname(idnt,gen_form(1),RUN_NUM,'$SET')
      call hbname(idnt,mctrue_form(1),MCELE,'$SET')
      call hbname(idnt,BGDTUP_form(1),EVTwant,'$SET')
    
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
            call bgdlumi97(histofill)

            WT_VTX = 1.0

	      if (histofill.eq.3) then
C ===========================================
C === GET A WEIGHT FOR VERTEX REWEIGHTING ===
C ===========================================
c               IF (YEAR.EQ.1996) THEN
c                 CALL VTX_WEIGHT(MC_ZV,WT_VTX,963,96,0)
c               call Rewt_Vtx_Isr(MC_ZV,0,WT_VTX,Ierr)
c               ELSEIF (YEAR.EQ.1997) THEN
                 CALL VTX_WEIGHT(MC_ZV,WT_VTX,971,97,0)
c               ENDIF
c ------------------------------------------
c Do lumi correction for ISR MC
c ------------------------------------------            
               call corrmclumig97         
		  end if

c ------------------------------------------
c General corrections and calculations F2ISR
c ------------------------------------------
            call corrections97
         	  electron_en = corrected_en
            electron_th = best_th

            call isr_calculate97

c ------------------------------------------
c Now apply cuts
c ------------------------------------------               
            acceptisr = 0
	      acceptempz = 0
		  calnorm = 0
		  kpaccept = 0

      call isr_cuts97(acceptisr, acceptempz,kpaccept,
     &                         calnorm,histofill)
c ------------------------------------------
c Calculate Weights
c ------------------------------------------           
		  if ((acceptisr.eq.1).or.(acceptempz.eq.1)) then

c ------------------------------------------
c Montecarlo
c ------------------------------------------
               if (histofill.eq.3) then
c	            if (year.eq.1996) then
c                    call Q2WTNEW96(q2_tru, mcwtq2, mcinput)
c	            else if (year.eq.1997) then
	               call Q2WTNEW97(q2_tru, mcwtq2, mcinput)
c	            end if
		        call f2wt97
	            bpcor = beamcorq(Zgamma,q2_el,2,ierr)
c			    bpcor = 1.0
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
            if (acceptisr.eq.1) then 	
	        call count_in_bins(histofill)
              call isr_histos97(histofill)
              count = count + 1 	    
	      end if

c ------------------------------------------
c Fill KP Histos
c ------------------------------------------
c      if (kpaccept.eq.1) then

c Calculate weights
c      if (acceptance.ne.0) then
c	   datawt = 1. / acceptance
c	else 
c	   datawt = 0.
c	end if
c	mcwt = mcwtq2 * f2weight * wt_vtx * bpcor

c Fill Histos		
c	temp = electron_en + elumig
c	if (histofill.eq.1) then
c	if ((electron_en.ge.4.).and.(electron_en.lt.7)) then
c	call hf1(500,temp,datawt)
c	else if((electron_en.ge.7.).and.(electron_en.lt.10)) then
c	call hf1(501,temp,datawt)
c	else if((electron_en.ge.10.).and.(electron_en.lt.13)) then
c	call hf1(502,temp,datawt)
c	else if((electron_en.ge.13.).and.(electron_en.lt.16)) then
c	call hf1(503,temp,datawt)
c	else if((electron_en.ge.16.).and.(electron_en.lt.19)) then
c	call hf1(504,temp,datawt)
c	else if((electron_en.ge.19.).and.(electron_en.lt.22)) then
c	call hf1(505,temp,datawt)
c	end if

c	else if (histofill.eq.3) then
c	if ((electron_en.ge.4.).and.(electron_en.lt.7)) then
c	call hf1(506,temp,mcwt)
c	else if((electron_en.ge.7.).and.(electron_en.lt.10)) then
c	call hf1(507,temp,mcwt)
c	else if((electron_en.ge.10.).and.(electron_en.lt.13)) then
c	call hf1(508,temp,mcwt)
c	else if((electron_en.ge.13.).and.(electron_en.lt.16)) then
c	call hf1(509,temp,mcwt)
c	else if((electron_en.ge.16.).and.(electron_en.lt.19)) then
c	call hf1(510,temp,mcwt)
c	else if((electron_en.ge.19.).and.(electron_en.lt.22)) then
c	call hf1(511,temp,mcwt)
c	end if

c	end if
c	end if

c ------------------------------------------
c Fill E-Pz Histos
c ------------------------------------------
		  if (acceptempz.eq.1) then
	         call empzfill97(histofill)
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
