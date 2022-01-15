c     ====================================================
       subroutine non_unif_cor(x,y,z,montecarlo,unif_cor)
c     ====================================================
c
c     non_unif_cor returns a correction factor for the measured electron
c     energy given the calorimeter x, y and z position of the electron.
c  
c     Outside of the region close to the cracks the correction value 
c     returned is 1.0, otherwise it varies between 1.0 and ~1.15).
c
c     For correct use of this non-uniformity correction routine it only
c     should be used together with RCALCORR version 6.0 or higher.
c
c     The correct calling sequence is:
c
c  in your eaze job: 
c 
c       1) RCALCORR version 6.0 (or higher)
c  
c  in the analysis job:
c
c       2) dead_material_energy
c       3) non_unif_cor
c     
c     The parameters needed for the correction are in the file 
c     NON_UNIF_PARA_2_0.INC which therefore has to be included.
c     Different parameters are used for data and Monte-Carlo.
c
c     In case of questions please send email to: deffner@desy.de 
c
c INPUT:
c
c     real x  = RCAL x position
c
c     real y  = RCAL y position
c
c     real z  = RCAL z position (only used to make sure no correction is
c               applied in the BCAL or FCAL region)
c     
c     logical montecarlo = true if Monte-Carlo     
c                          false if data
c
c OUTPUT:
c
c     real unif_cor = correction factor to be applied to measured electron
c                     energy  
c
c
C Author: Rolf Deffner  25-November-1998   Version 2.0
c
c  changes with respect to version 1.0:
c          
c   * always CAL electron position to be used (no more best_xyz needed)
c   * routine now factorizes with RCALCORR version 5.0 (and higher)   
c   * internally now same corr. factors for left and right CAL half
c
c
       implicit none
#include "non_unif_para_v2_0.inc"

       integer i, x_offset, y_offset
       real    x, y, z, unif_cor,xloss, yloss, xpar(5), ypar(5)
       logical montecarlo,apply_xcor, apply_ycor
       real    line, gauss, x_exp, y_exp        

c
c check whether z is in RCAL region
c
       if(z.gt.-140)then
        unif_cor = 1.
        return
       endif
c
c reset everything
c
      do i = 1, 5
       xpar(i) = 0. 
       ypar(i) = 0.
      enddo
      line = 0.
      gauss = 0.
      xloss = 1.
      yloss = 1.
      apply_xcor = .false.
      apply_ycor = .false.       
      x_offset = 999
      y_offset = 999
c
c =====
c DATA
c =====
c
      if(.NOT.montecarlo)then
c
c -------------------
c determine x offset (1 - 20)
c -------------------
c
        if(x.gt.-11.1.and.x.lt.-9.1.and.y.gt.0.0)  x_offset =   0
        if(x.gt.  9.1.and.x.lt.11.1.and.y.gt.0.0)  x_offset =   5
        if(x.gt.-11.1.and.x.lt.-9.1.and.y.lt.0.0)  x_offset =  10
        if(x.gt.  9.1.and.x.lt.11.1.and.y.lt.0.0)  x_offset =  15
c
c ------------------
c determine y offset
c ------------------
c
c both CAL halves (21 - 65)
c
       if(abs(x).gt.10.165)then
        if(y.gt.-71.5.and.y.lt.-65.5)       y_offset =  20
        if(y.gt.-52.0.and.y.lt.-46.0)       y_offset =  25
        if(y.gt.-32.5.and.y.lt.-26.5)       y_offset =  30
        if(y.gt.-12.0.and.y.lt. -8.0)       y_offset =  35
        if(y.gt. -2.3.and.y.lt.  1.7)       y_offset =  40
        if(y.gt.  7.5.and.y.lt. 11.5)       y_offset =  45
        if(y.gt. 26.1.and.y.lt. 32.1)       y_offset =  50
        if(y.gt. 45.7.and.y.lt. 51.7)       y_offset =  55
        if(y.gt. 65.5.and.y.lt. 71.5)       y_offset =  60
c
c module 12 (106 - 145)
c
       elseif(abs(x).lt.10.165)then
        if(y.gt.-66.5.and.y.lt.-60.5)       y_offset =  65  
        if(y.gt.-46.1.and.y.lt.-40.1)       y_offset =  70
        if(y.gt.-26.7.and.y.lt.-20.7)       y_offset =  75
        if(y.gt.-15.2.and.y.lt.-11.2)       y_offset =  80
        if(y.gt. 11.1.and.y.lt. 15.1)       y_offset =  85
        if(y.gt. 20.7.and.y.lt. 26.7)       y_offset =  90
        if(y.gt. 40.6.and.y.lt. 46.6)       y_offset =  95
        if(y.gt. 59.8.and.y.lt. 65.8)       y_offset = 100
       endif
c
c ===
c MC
c ===
c
c ------------------
c determine x offset (1 - 40)
c ------------------
c
      elseif(montecarlo)then
       if(y.gt.0.0)then
        if(x.gt.-11.4.and.x.lt. -8.4)  x_offset =   0
        if(x.gt.  8.4.and.x.lt. 11.4)  x_offset =   5
        if(x.gt.-30.8.and.x.lt.-27.8)  x_offset =  10
        if(x.gt. 27.8.and.x.lt. 30.8)  x_offset =  15
       elseif(y.lt.0.0)then 
        if(x.gt.-11.4.and.x.lt. -8.4)  x_offset =  20
        if(x.gt.  8.4.and.x.lt. 11.4)  x_offset =  25
        if(x.gt.-30.8.and.x.lt.-27.8)  x_offset =  30
        if(x.gt. 27.8.and.x.lt. 30.8)  x_offset =  35
       endif     
c
c ==================
c determine y offset
c ==================
c
c both CAL halves (41 - 85)
c
       if(abs(x).gt.10.165)then
        if(y.gt.-50.9.and.y.lt.-46.9)       y_offset =  40
        if(y.gt.-31.4.and.y.lt.-27.4)       y_offset =  45
        if(y.gt.-21.8.and.y.lt.-17.8)       y_offset =  50
        if(y.gt.-12.0.and.y.lt.- 8.0)       y_offset =  55
        if(y.gt. -2.3.and.y.lt.  1.7)       y_offset =  60
        if(y.gt.  7.4.and.y.lt. 11.4)       y_offset =  65
        if(y.gt. 17.1.and.y.lt. 21.1)       y_offset =  70
        if(y.gt. 26.8.and.y.lt. 30.8)       y_offset =  75
        if(y.gt. 46.6.and.y.lt. 50.6)       y_offset =  80
c
c module 12 (86 - 115)
c
       elseif(abs(x).lt.10.165)then
        if(y.gt.-65.0.and.y.lt.-59.0)       y_offset =  85  
        if(y.gt.-45.6.and.y.lt.-39.6)       y_offset =  90
        if(y.gt.-26.7.and.y.lt.-20.7)       y_offset =  95
        if(y.gt. 20.5.and.y.lt. 26.5)       y_offset = 100
        if(y.gt. 39.4.and.y.lt. 45.4)       y_offset = 105
        if(y.gt. 59.0.and.y.lt. 65.0)       y_offset = 110
       endif
      endif
c
c -----------------------------------------------------
c get the right parameters for the correction function
c -----------------------------------------------------
c
      if(.NOT.montecarlo)then
       if(x_offset.ne.999)then
        apply_xcor = .true. 
        do i = 1, 5
         xpar(i) =  DATA_EFRAC_COR(x_offset+i)
        enddo
       endif
       if(y_offset.ne.999)then
        apply_ycor = .true. 
        do i = 1, 5
         ypar(i) =  DATA_EFRAC_COR(y_offset+i)
        enddo
       endif
      elseif(montecarlo)then
       if(x_offset.ne.999)then
        apply_xcor = .true. 
        do i = 1, 5
         xpar(i) =  MC_EFRAC_COR(x_offset+i)
        enddo
       endif
       if(y_offset.ne.999)then
        apply_ycor = .true. 
        do i = 1, 5
         ypar(i) =  MC_EFRAC_COR(y_offset+i)
        enddo
       endif
      endif
c
c -----------------------------------------------------------------------
c     calculate x-dependent correction factor as the ratio of the fitted
c     line and the sum of the fitted gaussian and the line:
c          
c           xloss = line / (line + gauss) 
c -----------------------------------------------------------------------
c
       if(apply_xcor)then
        line = xpar(4)+ xpar(5)*x
        x_exp = -0.5*((x-xpar(2))/xpar(3))**2
        if(abs(x_exp).lt.10)then
         gauss = xpar(1)*exp(x_exp)
        else
         gauss = 0.
        endif
        xloss = line / (gauss + line)
       else
        xloss = 1.
       endif
c
c ----------------------------------------
c calculate y-dependent correction factor 
c ----------------------------------------
c
      if(apply_ycor)then
        line = ypar(4)+ ypar(5)*y
        y_exp = -0.5*((y-ypar(2))/ypar(3))**2
        if(abs(y_exp).lt.10)then
         gauss = ypar(1)*exp(y_exp)
        else
         gauss = 0.
        endif
        yloss = line / (gauss + line)
       else
        yloss = 1.
       endif 

      unif_cor = xloss*yloss 

      return
      end



