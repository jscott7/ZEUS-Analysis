	SUBROUTINE CELSWP
	implicit none
#include "partap.inc"
#include "caltru.inc"
#include "fmckin.inc"
#include "zdskey.inc"
C
C == CELSWP == swaps 3 pairs of BCAL EMC cells which were wrongly cabled.
C   Masahiro Kuze, 20 April 1997
C   Update: 1 July 1997, MK
C     Hardware fixed on 25 March.  This routine FIXES the DATA for the runs
C     before then.  It does nothing to the later runs, nor to MC.
C     If you like to do vice versa, change the lines between "C970701"
C ===THIS ROUTINE IS PRIVATE, NOT AN OFFICIAL STATEMENT FROM CAL GROUP ===
C
C  The problem exists at least from the beginning of 1993 running,
C  and was fixed in the access on 25 Mar 1997.
C  This routine fixes the mistake if applied to data event
C  (assuming it is not done in [re]processing, which is the case as of today)
C  and simulates the mistake if applied to MC event.
C
C  Input: none
C  Usage: CALL CELSWP before you look at CALTRU table.
C  Effect: It modifies contents of CALTRU (Energy, Imbalance and Time.)
C          A new cell entry may be created if it did not exist.  In this case
C          CConSa,CuPaOb link is INULL.  Recreate CConSa etc. yourself.
C  Caveat: If a cell was not in CALTRU, E=Imb=t(1)=t(2)=0.0 is assumed.
C          One has to do this game in RAW tape to be perfect.
C  CAUTION: This fix is not correct if any of the cells had imbalance=0.0
C   i.e. marked as a cell with a bad channel.  In this case the routine
C   modifies the table incorrectly and prints error message.  My guess is
C   this should not happen too often for 94-97 data.
C   For the 93 data, this does happen for one cell most of the time.
C   If you like to do it correctly for 93 data, please contact me.
C
C
	integer i, j, ii, jj, indcel, found(2)
	real el(2), eh(2), tl(2), th(2)
C handles (so far:-) 3 miscabled pairs.
	integer NPAIR
	parameter (NPAIR=3)
C list of the cells (Side information from Mike Jing)
	integer icel(2, NPAIR), side(NPAIR)
	data icel
C        cell #   ! Module Tower Cell Side-swapped(1=high, 0=low)
     +  /17458,   !   3      4    1     1 
     +   17462,   !   3      4    3     1
     +   25718,   !  19      8    3     0
     +   25720,   !  19      8    4     0
     +   32360,   !  32      7    4     1
     +   32370/   !  32      8    1     1
	data side /1,0,1/

C970701 from here

C If MC, return.

        if(coutab(FMCKIN).gt.0)return

C Last ep run with the mistake = 25497.  First run after fix = 25515.

        if(ZDSKEY_Nr1.gt.25500)return

C970701 ends here

	indcel=getind(CALTRU, 'Cellnr')

C loop over the NPAIR pairs.
	Do 111 i = 1, NPAIR
C get entries from CALTRU
	 Do 112, j = 1, 2
	  caltru_cellnr = icel(j, i)
	  call seltab(CALTRU, indcel, ii, jj)

	  if(jj.lt.ii)then
c entry not found.  Set all variables to 0.
	   el(j) = 0.
	   eh(j) = 0.
	   tl(j) = 0.
	   th(j) = 0.
	   found(j) = -1
	  else
c entry found - take 1st one.  If more entry, print error message.
	   if(jj.gt.ii)print *,'CELSWP: FATAL: double entry in CALTRU'
	   call fettab(CALTRU, indcel, ii)
	   el(j) = (caltru_e - caltru_imbal ) / 2.
	   eh(j) = (caltru_e + caltru_imbal ) / 2.
	   tl(j) = caltru_t(1)
	   th(j) = caltru_t(2)
	   found(j) = caltru_ID
c error message in case of bad cell.
c Please contact Masahiro.Kuze@desy.de if you see this message.
	   if(caltru_imbal.eq.0.0)print *, 
     +               'CELSWP: FATAL: bad cell. Result not correct'
c        Call PriTab(CALTRU,ID,MINC,MAXC,ALLCOL)
	  endif

112	 Continue

c do nothing if none of the 2 cells were present in caltru
	 if(found(1).eq.-1.and.found(2).eq.-1)goto 111

	 Do 113, j = 1, 2

c swap back, depending on the side which was wrongly cabled.
	  if(side(i).eq.0)then

	   caltru_e   = el(3-j) + eh(j)
	   caltru_imbal = eh(j) - el(3-j)
	   caltru_t(1) = tl(3-j)
	   caltru_t(2) = th(j)

	  else

	   caltru_e   = el(j) + eh(3-j)
	   caltru_imbal = eh(3-j) - el(j)
	   caltru_t(1) = tl(j)
	   caltru_t(2) = th(3-j)

	  endif

c replace entry, or newly insert 
c NB: It seems that, for an existing row, INULL does not overwrite.

	   caltru_cconsa = INULL
	   caltru_cupaob = INULL
	   caltru_cellnr = icel(j, i)

	   if(found(j).eq.-1)then
	    caltru_ID = NEXT
	    CALL INSENT(CALTRU)
	   else
	    caltru_ID = found(j)
	    CALL REPENT(CALTRU)
	   endif

113	 Continue

c        Call PriTab(CALTRU,ID,MINC,MAXC,ALLCOL)
111	Continue

	return
	end

