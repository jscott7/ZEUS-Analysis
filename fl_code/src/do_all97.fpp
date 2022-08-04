c ------------------------------
      subroutine do_all97(num)
c ------------------------------

      implicit none

#include "readdata.inc"
#include "local.inc"

      integer I,num
      integer choice, loop97
      real xmin,xmax,ymean
      loop97 = 2
	
c -----------------------------------------
c Read in Data
c -----------------------------------------
      call steer(loop97,num)
      write(*,*) '******************************'
      write(*,*) '***  97 Analysis selected  ***'	
      write(*,*) '***     Loop over Data     ***'
      write(*,*) '******************************'
      write(999,*) '******************************'
      write(999,*) '***  97 Analysis selected  ***'	
      write(999,*) '***     Loop over Data     ***'
      write(999,*) '******************************'

      datanorm = 0.

c Define bin
      fracdif = diffractiv
      q2min = 1.
      q2max = 30.
      ymin = ybinlow
      ymax = ybinhigh
      xmin = q2min / (90200. * ymin)
      xmax = q2max / (90200. * ymax)

      write(*,*) '********************************'
      write(*,*) '*** Q^2 and y bin boundaries ***'
      write(*,*) '********************************'
      write(999,*) '********************************'
      write(999,*) '*** Q^2 and y bin boundaries ***'
      write(999,*) '********************************'
 
      ymean = (10.**((log10(ymin) + log10(ymax))/2.))
      q2mean = (10.**((log10(q2min) + log10(q2max))/2.))
      xmean = q2mean/(ymean*90200.)  

      logxmean = log10(xmean)

      write(*,*) q2min, q2max, ymin, ymax
      write(999,*) q2min, q2max, ymin, ymax

c -----------------------------------------
c 1997 Data
c ------------------------------------------
      call readdat(45)
      do I = 1, NR_INFILES
         call ANEVENT97(FILEN(I)(1:LENOCC(FILEN(I))),10,1)
      end do
      close (45)
c -------------------------------------------------
c Read in Background x1 - calculate normalisation
c -------------------------------------------------
      write(*,*) '******************************'
      write(*,*) '*** Background Loop #1     ***'
      write(*,*) '*** Calculate Normalisation***'
      write(*,*) '******************************'
      write(999,*) '******************************'
      write(999,*) '*** Background Loop #1     ***'
      write(999,*) '*** Calculate Normalisation***'
      write(999,*) '******************************'

      backnorm = 0.
      baknorm = 0.
      call readdat(46)
      do I = 1, NR_INFILES
         call ANEVENT97(FILEN(I)(1:LENOCC(FILEN(I))),20,20)
      end do
      close(46)

      baknorm = datanorm / backnorm
      write(*,*) 'Normalisation Factor', baknorm
      write(999,*) 'Normalisation Factor', baknorm
      write(*,*) '# background events including E-Pz > 62', backnorm
      write(999,*) '# background events including E-Pz > 62', backnorm

c -------------------------------------------------
c Read in Background x2 - Fill histos
c -------------------------------------------------
      write(*,*) '******************************'
      write(*,*) '*** Background loop #2     ***'
      write(*,*) '*** Fill Histograms        ***'
      write(*,*) '******************************' 
      write(999,*) '******************************'
      write(999,*) '*** Background loop #2     ***'
      write(999,*) '*** Fill Histograms        ***'
      write(999,*) '******************************' 

      call readdat(46)
      do I = 1, NR_INFILES
         call ANEVENT97(FILEN(I)(1:LENOCC(FILEN(I))),20,2)
      end do
      close(46)

c -----------------------------------------
c Read in Montecarlo
c -----------------------------------------
      write(*,*) '******************************'
      write(*,*) '***     Loop over MC       ***'
      write(*,*) '******************************'
      write(999,*) '******************************'
      write(999,*) '***     Loop over MC       ***'
      write(999,*) '******************************'

      mcinput = 0
      call readdat(47)
      do I = 1, NR_INFILES
         mcinput = I
         call ANEVENT97(FILEN(I)(1:LENOCC(FILEN(I))),10,3)
      end do
      close(47)  
c ------------------------------------------
c Now do F2 calculation
c ------------------------------------------
      call calculate_f2
    
      return
      end

