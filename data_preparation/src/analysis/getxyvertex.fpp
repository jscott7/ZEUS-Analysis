      subroutine GetXyVertex( Run, X, xerr, Y, yerr ,ierr )
C     -------------------------------------------------
C     Return X and Y of average vertex given run number
C     -------------------------------------------------
C
C     Input: 
C       integer run  : for MC enter run=1 and I will check if
C                      you are using 96 or 97 MC from the tape name.
C
C     Output:
C         real X,Zerr    - (out) Average X vertex and error
C         real Y,Yerr    - (out) Average Y vertex and error
C         integer ierr
C                        0 - All is okay
C                        1 - Run out of range 
C                        2 - MC file name not understood
C
C  Version 1.0 Mike Wodarczyk Feb 9,1999
C
      IMPLICIT NONE
#include "partap.inc"
#include "zrinpt.inc"

      Integer Run
      Real X, Y
      real xerr,yerr
      integer ierr

      integer noevt
      logical firsterr
      data firsterr /.TRUE. /
      save firsterr

      ierr = 0
      if ( run .gt. 1 ) then
         CALL  XYVTX(Run,X,XERR,Y,YERR,NOEVT)
         if ( noevt .ge. 0 ) then
            ierr = 0
         else
            if ( firsterr ) then
               write(*,*) 'ERROR: getxyvertex doesnt know run  '
               write(*,*) 'run = ',run
               write(*,*) 'Using 1997 MC xy vertex'
               write(*,*) 'getxyvertex only complains once so '//
     &              'listen to me.'
               firsterr = .FALSE.
            endif
            ierr = 1
            x = -0.08186        ! 1997 genvtx
            xerr = 0.0001
            y = -0.03996        ! 1997 genvtx
            yerr = 0.0001
         endif
      else
C MC
        if ( ZRINPT_NAME(1)(1:6).eq.'ZEUSMC' ) then
C...  Name is a normal MC file   ZEUSMC.A6WP613.
           if ( ZRINPT_NAME(1)(12:12).eq.'6' ) then
              x = -0.0333       ! 1996 genvtx
              xerr = 0.0001
              y = -0.0381       ! 1996 genvtx
              yerr = 0.0001
              ierr = 0
           else if ( ZRINPT_NAME(1)(12:12).eq.'7' ) then
              x = -0.08186      ! 1997 genvtx
              xerr = 0.0001
              y = -0.03996      ! 1997 genvtx
              yerr = 0.0001
              ierr = 0
           else
              if ( firsterr ) then
                 write(*,*) 'ERROR: getxyvertex doesnt understand '// 
     &                'MC version '
                 write(*,*) 'MC File = '//ZRINPT_NAME(1)
                 write(*,*) 'Using 1997 xy vertex'
                 write(*,*) 'getxyvertex only complains once so '//
     &                'listen to me.'
                 firsterr = .FALSE.
              endif
              x = -0.08186      ! 1997 genvtx
              xerr = 0.0001
              y = -0.03996      ! 1997 genvtx
              yerr = 0.0001
              ierr = 2
           endif
        endif
      endif
C
C

      RETURN
      END

C
      SUBROUTINE XYVTX(RUNNO,OXV,OXVERR,OYV,OYVERR,ONOEVT)
C====================================================================
C     this subroutine needs rvtx96.dat
C                      &    rvtx97.dat
C
C     INPUT:    RUNNO
C     OUTPUT:   OXV,OZVERR       XVTX AND ERROR
C               OYV,OZVERR       YVTX AND ERROR
C               ONOEVT           NUMBER OF EVENTS USED FOR THE
C                                VERTEX CALCULATION
C
C     if the run number does not exist, onoevt will be -99
C
C     OCR Feb99
C     version 1.0  will be changed soon      02/02/99
C
C====================================================================
      IMPLICIT NONE

      LOGICAL FIRST
      DATA    FIRST /.TRUE./
      INTEGER I, NOR
      INTEGER COMPRUN
      DATA    COMPRUN /0/
      SAVE    COMPRUN
C..  INPUT
      INTEGER RUNNO
C .. OUTPUT:
      REAL    OXV,OXVERR,OYV,OYVERR
      REAL    OXV_S,OXVERR_S,OYV_S,OYVERR_S
      DATA    OXV_S /0.0/, OXVERR_S /0.0/ ,OYV_S /0.0/, OYVERR_S /0./
      SAVE    OXV_S,OXVERR_S,OYV_S,OYVERR_S      

      INTEGER ONOEVT,ONOEVT_S
      DATA    ONOEVT_S /0/
      SAVE    ONOEVT_S

      REAL    RUN_NO(1000),NOEVT(1000)
      REAL    XV(1000),XVERR(1000)
      REAL    YV(1000),YVERR(1000)

      IF (RUNNO .EQ. COMPRUN) THEN
         OXV    = OXV_S
         OXVERR = OXVERR_S
         OYV    = OYV_S
         OYVERR = OYVERR_S
         ONOEVT = ONOEVT_S
      ELSE
         OXV    = -99.9
         OXVERR = -99.9
         OYV    = -99.9
         OYVERR = -99.9
         ONOEVT = -99
         COMPRUN = RUNNO
         IF (RUNNO .GT. 20000 .AND. RUNNO .LT. 23000) THEN
            OPEN (UNIT=21,FILE='rvtx96.dat',STATUS='OLD')
         ELSEIF (RUNNO .GT. 24999 .AND. RUNNO .LT. 28000) THEN
            OPEN (UNIT=21,FILE='rvtx97.dat',STATUS='OLD')
         ELSE
            WRITE(*,*) 'UNKNOWN RUN NUMBER'
            RETURN
         ENDIF
         NOR = 0
         DO I = 1,1000
            NOR = NOR + 1
            READ (21,*,ERR=33) RUN_NO(I),XV(I),XVERR(I),
     &                                   YV(I),YVERR(I),NOEVT(I)
            IF (RUNNO.EQ.RUN_NO(I)) THEN
               OXV    = XV(I)
               OXVERR = XVERR(I)
               OYV    = YV(I)
               OYVERR = YVERR(I)
               ONOEVT = INT(NOEVT(I))
               OXV_S    = XV(I)
               OXVERR_S = XVERR(I)
               OYV_S    = YV(I)
               OYVERR_S = YVERR(I)
               ONOEVT_S = INT(NOEVT(I))
            ENDIF
         ENDDO
 33      CONTINUE
         CLOSE(21)
         WRITE(*,*)'XYVTX: RUN_NO, X, Y', runno, oxv, oyv
      ENDIF

      RETURN
      END
