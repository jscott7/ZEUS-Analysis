c     =====================================================
       subroutine non_unif_corv1(x,y,z,montecarlo,unif_cor)
c     =====================================================
c
c     non_unif_cor returns a correction factor for the measured electron
c     energy given the x, y and z position of the electron.
c  
c     The correction functions were obtained by parametrizing the dips 
c     in the Emeas/Eda-distributions close to the RCAL-(module/tower/cell-)
c     cracks in x and y. 
c     When determing the correction factors the UNCORRECTED CAL ENERGY 
c     was used, i.e. NO SRTD- OR PRES-CORRECTION was applied.
c  
c     Outside of the region close to the cracks the correction value 
c     returned is 1.0, otherwise it varies between 1.0 and ~1.15).
c
c     Different components were used for a RCAL radius R > 50 cm 
c     (CTD track extrapolation used) and for R < 50 cm (SRTD/HES/RCAL).
c     So in order to use the correction routine correctly outside a RCAL 
c     radius of 50 cm the position of the electron track should be used 
c     as input parameter. For R < 50 cm the routine BEST_XYZ returns the
c     x, y and z coordinate of the electron using either SRTD, HES or CAL.
c     This routine is a version of Mike's best_theta routine, the only 
c     difference is that it returns x, y and z instead of theta.
c
c     The parameters needed for the correction are in the file 
c     NON_UNIF_PARA.INC which therefore has to be included.
c     Different parameters are used for data and Monte-Carlo.
c
c     In case of questions please send email to: deffner@desy.de 
c
c INPUT:
c
c     real x  = x at RCAL
c
c     real y  = y at RCAL 
c
c     real z  = z at RCAL  (only used to make sure no correction is applied
c                           in the BCAL or FCAL region)
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
C Author: Rolf Deffner  03-July-1998   Version 1.0
c
c

       implicit none
#include"non_unif_para.inc"

       integer i, x_offset, y_offset
       real    x, y, z, r, unif_cor,xloss, yloss, xpar(5), ypar(5)
       logical montecarlo,apply_xcor, apply_ycor
       real    line, gauss        

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
c find out in which region we are
c
      r = sqrt(x**2 + y**2) 
c
c DATA
c
      if(.NOT.montecarlo)then
c determine x offset
       if(r.gt.50.0)then
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.gt.0.0)  x_offset =   0
        if(x.gt.  8.0.and.x.lt.12.0.and.y.gt.0.0)  x_offset =   5
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.lt.0.0)  x_offset =  10
        if(x.gt.  8.0.and.x.lt.12.0.and.y.lt.0.0)  x_offset =  15
       elseif(r.lt.50.0)then
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.gt.0.0)  x_offset =  20
        if(x.gt.  8.0.and.x.lt.12.0.and.y.gt.0.0)  x_offset =  25
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.lt.0.0)  x_offset =  30
        if(x.gt.  8.0.and.x.lt.12.0.and.y.lt.0.0)  x_offset =  35
       endif  
c ==================
c calculate y offset
c ==================
c R > 50 cm and left CAL half
       if(r.gt.50.0.and.x.lt.-10.0)then
        if(y.gt.-75.0.and.y.lt.-60.0)       y_offset =  40
        if(y.gt.-55.0.and.y.lt.-40.0)       y_offset =  45
        if(y.gt.-35.0.and.y.lt.-23.0)       y_offset =  50
        if(y.gt.-23.0.and.y.lt.-15.0)       y_offset =  55
        if(y.gt.-15.0.and.y.lt. -5.0)       y_offset =  60
        if(y.gt. -5.0.and.y.lt.  5.0)       y_offset =  65
        if(y.gt.  5.0.and.y.lt. 15.0)       y_offset =  70
        if(y.gt. 15.0.and.y.lt. 22.0)       y_offset =  75
        if(y.gt. 25.0.and.y.lt. 35.0)       y_offset =  80
        if(y.gt. 40.0.and.y.lt. 55.0)       y_offset =  85
        if(y.gt. 55.0.and.y.lt. 75.0)       y_offset =  90
c R > 50 cm and right CAL half
       elseif(r.gt.50.0.and.x.gt.10.0)then
        if(y.gt.-75.0.and.y.lt.-60.0)       y_offset =  95
        if(y.gt.-55.0.and.y.lt.-40.0)       y_offset = 100
        if(y.gt.-35.0.and.y.lt.-23.0)       y_offset = 105
        if(y.gt.-23.0.and.y.lt.-15.0)       y_offset = 110
        if(y.gt.-15.0.and.y.lt. -5.0)       y_offset = 115
        if(y.gt. -5.0.and.y.lt.  5.0)       y_offset = 120
        if(y.gt.  5.0.and.y.lt. 15.0)       y_offset = 125
        if(y.gt. 15.0.and.y.lt. 25.0)       y_offset = 130
        if(y.gt. 25.0.and.y.lt. 35.0)       y_offset = 135
        if(y.gt. 40.0.and.y.lt. 55.0)       y_offset = 140
        if(y.gt. 55.0.and.y.lt. 75.0)       y_offset = 145
c R < 50 cm, left CAL half
       elseif(r.lt.50.0.and.x.lt.-10.0)then
        if(y.gt.-33.0.and.y.lt.-25.0)       y_offset = 150 
c no correction in this region for data     y_offset = 155
        if(y.gt.-17.0.and.y.lt. -2.0)       y_offset = 160
        if(y.gt. -2.0.and.y.lt.  2.0)       y_offset = 165
        if(y.gt.  0.0.and.y.lt. 18.0)       y_offset = 170
c no correction in this region for data     y_offset = 175
        if(y.gt. 23.0.and.y.lt. 33.0)       y_offset = 180
c R < 50 cm, right CAL half
       elseif(r.lt.50.0.and.x.gt.10.0)then
        if(y.gt.-35.0.and.y.lt.-25.0)       y_offset = 185
c no correction in this region for data     y_offset = 190 
        if(y.gt.-14.0.and.y.lt. -6.0)       y_offset = 195
        if(y.gt. -3.0.and.y.lt.  3.0)       y_offset = 200
        if(y.gt.  4.0.and.y.lt. 14.0)       y_offset = 205
c no correction in this region for data     y_offset = 210 
        if(y.gt. 25.0.and.y.lt. 32.0)       y_offset = 215
c module 12, R < 50 cm
       elseif(r.lt.50.0.and.abs(x).lt.10.0)then
        if(y.gt.-50.0.and.y.lt.-37.0)       y_offset = 220  
        if(y.gt.-27.0.and.y.lt.-20.0)       y_offset = 225
        if(y.gt.-16.0.and.y.lt.-10.0)       y_offset = 230
        if(y.gt. 10.0.and.y.lt. 16.0)       y_offset = 235
        if(y.gt. 21.0.and.y.lt. 26.0)       y_offset = 240
        if(y.gt. 35.0.and.y.lt. 50.0)       y_offset = 245
       endif
c
c MC
c
      elseif(montecarlo)then
       if(r.gt.50.0)then
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.gt.0.0)  x_offset =   0
        if(x.gt.  8.0.and.x.lt.12.0.and.y.gt.0.0)  x_offset =   5
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.lt.0.0)  x_offset =  10
        if(x.gt.  8.0.and.x.lt.12.0.and.y.lt.0.0)  x_offset =  15
       elseif(r.lt.50.0)then
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.gt.0.0)  x_offset =  20
        if(x.gt.  8.0.and.x.lt.12.0.and.y.gt.0.0)  x_offset =  25
        if(x.gt.-12.0.and.x.lt.-8.0.and.y.lt.0.0)  x_offset =  30
        if(x.gt.  8.0.and.x.lt.12.0.and.y.lt.0.0)  x_offset =  35
       endif  
c ==================
c calculate y offset
c ==================
c R > 50 cm and left CAL half
       if(r.gt.50.0.and.x.lt.-10.0)then
        if(y.gt.-75.0.and.y.lt.-60.0)       y_offset =  40
        if(y.gt.-55.0.and.y.lt.-40.0)       y_offset =  45
        if(y.gt.-35.0.and.y.lt.-23.0)       y_offset =  50
        if(y.gt.-23.0.and.y.lt.-15.0)       y_offset =  55
        if(y.gt.-15.0.and.y.lt. -5.0)       y_offset =  60
        if(y.gt. -5.0.and.y.lt.  5.0)       y_offset =  65
        if(y.gt.  5.0.and.y.lt. 15.0)       y_offset =  70
        if(y.gt. 15.0.and.y.lt. 25.0)       y_offset =  75
        if(y.gt. 25.0.and.y.lt. 35.0)       y_offset =  80
        if(y.gt. 40.0.and.y.lt. 55.0)       y_offset =  85
        if(y.gt. 55.0.and.y.lt. 75.0)       y_offset =  90
c R > 50 cm and right CAL half
       elseif(r.gt.50.0.and.x.gt.10.0)then
        if(y.gt.-75.0.and.y.lt.-60.0)       y_offset =  95
        if(y.gt.-55.0.and.y.lt.-40.0)       y_offset = 100
        if(y.gt.-35.0.and.y.lt.-23.0)       y_offset = 105
        if(y.gt.-23.0.and.y.lt.-15.0)       y_offset = 110
        if(y.gt.-15.0.and.y.lt. -5.0)       y_offset = 115
        if(y.gt. -5.0.and.y.lt.  5.0)       y_offset = 120
        if(y.gt.  5.0.and.y.lt. 15.0)       y_offset = 125
        if(y.gt. 15.0.and.y.lt. 25.0)       y_offset = 130
        if(y.gt. 25.0.and.y.lt. 35.0)       y_offset = 135
        if(y.gt. 40.0.and.y.lt. 55.0)       y_offset = 140
        if(y.gt. 55.0.and.y.lt. 75.0)       y_offset = 145
c R < 50 cm, left CAL half
       elseif(r.lt.50.0.and.x.lt.-10.0)then
        if(y.gt.-33.0.and.y.lt.-25.0)       y_offset = 150 
        if(y.gt.-25.0.and.y.lt.-15.0)       y_offset = 155
        if(y.gt.-14.0.and.y.lt. -6.0)       y_offset = 160
        if(y.gt. -2.0.and.y.lt.  2.0)       y_offset = 165
        if(y.gt.  6.0.and.y.lt. 12.0)       y_offset = 170
        if(y.gt. 16.0.and.y.lt. 23.0)       y_offset = 175
        if(y.gt. 24.0.and.y.lt. 32.0)       y_offset = 180
c R < 50 cm, right CAL half
       elseif(r.lt.50.0.and.x.gt.10.0)then
        if(y.gt.-32.0.and.y.lt.-26.0)       y_offset = 185
        if(y.gt.-25.0.and.y.lt.-15.0)       y_offset = 190
        if(y.gt.-14.0.and.y.lt. -6.0)       y_offset = 195
        if(y.gt. -2.0.and.y.lt.  2.0)       y_offset = 200
        if(y.gt.  6.0.and.y.lt. 12.0)       y_offset = 205
        if(y.gt. 15.0.and.y.lt. 24.0)       y_offset = 210
        if(y.gt. 24.0.and.y.lt. 32.0)       y_offset = 215
c module 12, R < 50 cm
       elseif(r.lt.50.0.and.abs(x).lt.10.0)then
        if(y.gt.-46.0.and.y.lt.-33.0)       y_offset = 220  
        if(y.gt.-28.0.and.y.lt.-20.0)       y_offset = 225
        if(y.gt.-16.0.and.y.lt.-10.0)       y_offset = 230
        if(y.gt. 12.0.and.y.lt. 16.0)       y_offset = 235
        if(y.gt. 20.0.and.y.lt. 27.0)       y_offset = 240
        if(y.gt. 37.0.and.y.lt. 48.0)       y_offset = 245
       endif
      endif
c
c get the right parameters for the correction function
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
c     calculate x-dependent correction factor as the ratio of the fitted
c     line and the sum of the fitted gaussian and the line:
c          
c           xloss = line / (line + gauss) 
c
       if(apply_xcor)then
        line = xpar(4)+ xpar(5)*x
        gauss = xpar(1)*exp(-0.5*((x-xpar(2))/xpar(3))**2)
        xloss = line / (gauss + line)
       else
        xloss = 1.
       endif 
c
c calculate y-dependent correction factor 
c
      if(apply_ycor)then
        line = ypar(4)+ ypar(5)*y
        gauss = ypar(1)*exp(-0.5*((y-ypar(2))/ypar(3))**2)
        yloss = line / (gauss + line)
       else
        yloss = 1.
       endif 

      unif_cor = xloss*yloss 

      return
      end


