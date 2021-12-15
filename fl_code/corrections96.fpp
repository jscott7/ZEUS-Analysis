c --------------------------------
      subroutine corrections96
c --------------------------------

      implicit none

#include "common96.inc"
#include "local.inc"
#include "constant.inc"

      integer i, j
      real unif_cor
      real dead_mat_en
      integer tvctpar, mode, dmemode
      real cal(4), srtd(4), hes(4), Tend(5)
      real vtx(3)
      real dum1, dum2
	real nonunif_en
	integer whatToUse, Isys

      whatToUse = 0
	Isys = 0

	dead_mat_en = 0.
	nonunif_en = 0.	
	unif_cor = 0.

	do i = 1,3 
	vtx(i) = 0.
	end do

      do j = 1, 4
         cal(j) = 0.
         srtd(j) = 0.
         hes(j) = 0.
	   Tend(j) = 0.
      end do

      Tend(5) = 0.

      
      cal(1) = cal_xp
      cal(2) = cal_yp
      cal(3) = cal_zp
      cal(4) = cal_ee

      srtd(1) = srtd_x
      srtd(2) = srtd_y
      srtd(3) = srtd_z
c -----------------------------------
c SRTD energy in MIP, sum of 2 layers
c -----------------------------------
      srtd(4) = srtd_e
      
c      hes(1) = hes_x
c      hes(2) = hes_y
c      hes(3) = -99.
c      hes(4) = hes_e
	     
c      Tend(1) = trk_x
c      Tend(2) = trk_y
c      Tend(3) = trk_z
c      Tend(4) = trk_p
c      Tend(5) = trk_d

	hes(1) = -9999.
	hes(2) = -9999.
	hes(3) = -9999.
	hes(4) = -9999.

	tend(1) = -9999.
	tend(2) = -9999.
	tend(3) = -9999.
	tend(4) = -9999.
	tend(5) = -9999.

c ---------------------------------------------
c Parameters for best theta :
c tvctpar not used -> set to -1
c mode = 1 (just calculate best theta and quit)
c ---------------------------------------------
      tvctpar = -1
      mode = 1
c ---------------------------------------------
c Parameters for dead material energy correction :
c dmemode = 0 (Apply nominal correction)
c ---------------------------------------------
      dmemode = 0
      
	if (VCT_XVC.gt.-300.) then
	vtx(1) = VCT_XVC
	else
	vtx(1) = 0.
	end if

	if (VCT_YVC.gt.-300.) then
	vtx(2) = VCT_YVC
	else
	vtx(2) = 0.
	end if
	 
	if (VCT_ZVC.gt.-300.) then
	vtx(3) = VCT_ZVC
	else
	vtx(3) = 0.
	end if


c	call bestvtx(VCT_ZVC, NVTRKC, CHVCC,TRK_V,FCAL_VTX,FCAL_VTXE,
c     & cal(1),cal(2),cal(3),zgamma,vtx(3),whatToUse)

      call best_xyz(bestx, besty, bestz, cal,
     &     srtd, hes, Tend, vtx, montecarlo)
c -------------------------------------------------
c New version - but not final values.
c -------------------------------------------------
c      call best_position(best_th, bestx, besty, bestz,
c     & cal, srtd, hes, Tend, vtx, montecarlo, run_num, whatToUse,
c     & Isys)
		
      call best_theta(best_th, cal, srtd, hes, trk_t, tvctpar,
     &     vtx, mode, montecarlo) 

c ------------------------------------------------------
c Old non unif correction - using best XYZ
c COMMENT OUT THE FOLLOWING LINES WHEN USING NEW VERSION
c ------------------------------------------------------
c      call non_unif_corv1(bestx, besty, bestz, montecarlo, unif_cor)
      
c 	nonunif_en = cal_ee * unif_cor

c      call dead_material_energy(corrected_en,nonunif_en,cal,srtd, 
c     &     prse2, montecarlo, dmemode)

c ------------------------------------------------------
c New Non unif correction - don't need best XYZ
c ** Check later with old version **
c should be used with RCALCOR version 6.0 or above
c COMMENT OUT THE FOLLOWING LINES WHEN USING OLD VERSION
c ------------------------------------------------------
	call dead_material_energy(dead_mat_en,cal_ee,cal,srtd, 
     &     prse2, montecarlo, dmemode,srtdpressc)

      call non_unif_cor(cal_xp, cal_yp, cal_zp,
     &        montecarlo, unif_cor)
  
      corrected_en = dead_mat_en * unif_cor
      
      electron_en = corrected_en * elenscale

      electron_th = best_th

      return
      end








