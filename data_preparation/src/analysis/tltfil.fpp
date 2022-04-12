
C     ==================
      SUBROUTINE TLTBITS
C     ==================

       IMPLICIT NONE
c -------
c Author : Gareth Howell 29/5/97
c Routine to transfer TLT subtriggers into local logical arrays.
c -------
c
#include "partap.inc"
#include "zrevt.inc"
#include "tltevt.inc"
#include "common.inc"

       INTEGER i,JBIT, iT
c
       call vzero(TLT,15)
       if(coutab(TLTEVT).gt.0)then
          call FETTAB(TLTEVT,ID,1)
          do iT = 1, 15
             TLT(IT) = TLTEVT_Subtrg(IT)
          end do
       end if
             
c        IF (COUTAB(TLTEVT).GT.0) THEN
c          CALL FETTAB(TLTEVT,ID,1)
c          DO i=1,32
c            TLTSUB1(i)   = JBIT(TLTEVT_Subtrg(1),i).eq.1
c            TLTSUB3(i)   = JBIT(TLTEVT_Subtrg(3),i).eq.1
c            TLTSUB4(i)   = JBIT(TLTEVT_Subtrg(4),i).eq.1
c            TLTSUB5(i)   = JBIT(TLTEVT_Subtrg(5),i).eq.1
c            TLTSUB6(i)   = JBIT(TLTEVT_Subtrg(6),i).eq.1
c            TLTSUB7(i)   = JBIT(TLTEVT_Subtrg(7),i).eq.1
c            TLTSUB8(i)   = JBIT(TLTEVT_Subtrg(8),i).eq.1
c            TLTSUB9(i)   = JBIT(TLTEVT_Subtrg(9),i).eq.1
c            TLTSUB10(i)  = JBIT(TLTEVT_Subtrg(10),i).eq.1
c            TLTSUB11(i)  = JBIT(TLTEVT_Subtrg(11),i).eq.1
c            TLTSUB12(i)  = JBIT(TLTEVT_Subtrg(12),i).eq.1
c            TLTSUB13(i)  = JBIT(TLTEVT_Subtrg(13),i).eq.1
c          ENDDO
c        ELSE
c          write(6,*) ' TLTBITS error : TLTEVT empty RUN/EVT'
c     +    ,zrevt_runnr,zrevt_evtnr(3)
c
c         DO i=1,32
c            TLTSUB1(i)   = .false.
c            TLTSUB3(i)   = .false.
c            TLTSUB4(i)   = .false.
c            TLTSUB5(i)   = .false.
c            TLTSUB6(i)   = .false.
c            TLTSUB7(i)   = .false.
c            TLTSUB8(i)   = .false.
c            TLTSUB9(i)   = .false.
c            TLTSUB10(i)  = .false.
c            TLTSUB11(i)  = .false.
c            TLTSUB12(i)  = .false.
c            TLTSUB13(i)  = .false.
c          ENDDO
c        ENDIF
       END

