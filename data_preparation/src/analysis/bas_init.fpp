C     ===================
      SUBROUTINE BAS_INIT(MC)
C     ===================
C
C --- Basic Initialization ....
C
C     A.Quadt, May 1997.
C     -----------------------------------

      IMPLICIT NONE

#include "partap.inc"
c#include "ctrlcm.inc"
#include "fmckin.inc"
#include "zrevt.inc"
#include "common.inc"

      INTEGER      IERR
      character*4 genkey
      logical MC

C ----------------------------------------------------------------
C --- Load CAL geometry for use with electron finder routines. ---
C ----------------------------------------------------------------
      CALL CCGEOM(IERR)

      WRITE(*,*) ' '
      WRITE(*,*) ' '

      IF (IERR.EQ.0) THEN
          WRITE(*,*) ' ----------------------------------------- '
          WRITE(*,*) ' ---- CAL GEOMETRY SUCCESSFULLY LOADED --- '
          WRITE(*,*) ' ----------------------------------------- '
      ELSEIF (IERR.NE.0) THEN
          WRITE(*,*) ' ----------------------------------- '
          WRITE(*,*) ' ---- CAL GEOMETRY LOAD ERROR --- ',IERR
          WRITE(*,*) ' ----------------------------------- '
      ENDIF
      
      WRITE(*,*) ' '
      WRITE(*,*) ' '
 
c Load SRTD geometry setup
      call srgeom(ierr)
      if (ierr.ne.0) then
         write(*,*) 'SRTD geometry loading error:',ierr
      end if

      if ((zrevt_time(1).gt.19950000).and.
     &     (zrevt_time(1).lt.19960000)) then
         year = 1995
      elseif ((zrevt_time(1).gt.19960000).and.
     &        (zrevt_time(1).lt.19970000)) then
         year = 1996
      elseif ((zrevt_time(1).gt.19970000).and.
     &        (zrevt_time(1).lt.19980000)) then
         year = 1997
      elseif ((zrevt_time(1).gt.19980000).and.
     &        (zrevt_time(1).lt.19990000)) then
         year = 1998
      else
         write(*,*)'Could not set year, using 1995'
         year = 1995
      end if

      call defgen(genkey)
      if (genkey .eq. 'DATA') then
        MC = .false.
      else 
	MC = .true.
      end if   
      
      IF (.not.MC) THEN
        WRITE(*,*) ' '
        WRITE(*,*) ' ======== REAL DATA ANALYSIS SELECTED ======== '
        WRITE(*,*) ' '
 
      ELSEIF (MC) THEN
        WRITE(*,*) ' '
        WRITE(*,*) ' ======== MC ANALYSIS SELECTED ======== '
        WRITE(*,*) ' '
      ENDIF

      CALL ZESMAF        !!! Get B-Field
 
      IF (MC) THEN
          WRITE(*,*) ' ---------- INITIALISING O1RECON -------- '
          CALL O1INIT
          CALL O1ZGID
      ENDIF

      RETURN
      END








