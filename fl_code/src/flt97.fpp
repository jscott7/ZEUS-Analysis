C
C     ---------------------------------
      SUBROUTINE TRIGGERCUT97(TFLAG,IERR)
C     ---------------------------------
      implicit none

#include "local.inc"
#include "common97.inc"

      REAL remcthval,isoeval,remcval,caleval

      LOGICAL SRTD47_t_vgd,SRTD56_t_vgd,SRTD47_t_bad,SRTD56_t_bad,
     & SRTD_HIT,REMCTH,REMC,ISOE,SRTDGOOD,CALE

c      INCLUDE 'f2isr.inc'
      integer ierr
      logical trigfirst,tflag
      save trigfirst
      data trigfirst /.true./

      remcval=2032
      caleval=464
      remcthval=3750
      isoeval=0

      if (trigfirst) then
         WRITE (*,*) ' ****************************'
         WRITE (*,*) ' '
         WRITE (*,*) '     TRIGGER CUTS  '
         WRITE (*,*) ' '
         WRITE (*,*) ' ****************************'
         WRITE (999,*) ' ****************************'
         WRITE (999,*) ' '
         WRITE (999,*) '     TRIGGER CUTS  '
         WRITE (999,*) ' '
         WRITE (999,*) ' ****************************'

         trigfirst=.false.
      endif
C
C    RCAL_isoe*RE_th_gSRTD2*96g(2032,3750,464)
C
C    RCAL_isoe*RE_th_gSRTD2*96g
C
C    id: 100 
C    logic: Reg_RCAL_isoe &&
C    _____ ( RCAL_EMC_E >= $1 || REMCth >= $2 ||
C    _______ ( CAL_E3 >= $3 && SRTD_good ) ) 
C    veto: C5v, VWiv, VWov, SRTD95v2 
C
      SRTD_HIT=.false.
      REMCTH=.false.
      REMC=.false.
      ISOE=.false.
      SRTDGOOD=.false.
      CALE=.false.
C
      if (trigdat(10).gt.remcval) then
       remc=.true.
      endif
C
      if (trigdat(11).gt.caleval) then
       cale=.true.
      endif
C
      if (trigdat(12).gt.remcthval) then
       remcth=.true.
      endif
C
      if (trigdat(15).gt.isoeval) then
       isoe=.true.
      endif
C
c      if (srtd_mip(1).gt.0.0) then
	if (srtd_e.gt.0.0) then
       srtdgood=.true.
      endif
C
      if (isoe.and.(remcth.or.remc.or.(cale.and.srtdgood))) then
       tflag=.true.
      endif

      RETURN
      END
C



