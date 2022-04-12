C   ------------------------------------------------
      SUBROUTINE sprkrej
C   ------------------------------------------------
c   This routine loops over the CALTRU table and
c   removes entries for cell 11052 from the CALTRU 
c   table if the imbalance (caltru_imbal/caltru_e) 
c   is greater than 0.8
c
C    OCR 01/07/98


      IMPLICIT NONE

#include "partap.inc"
#include "caltru.inc"
#include "fmckin.inc"
#include "zdskey.inc"

      INTEGER   NRCELLS, I
      LOGICAL   PRINTTAB, FIRST
      DATA      FIRST /.TRUE./

      IF (FIRST) THEN
        FIRST = .FALSE.
        WRITE(*,*) '======== SPRKREJ ========='
      ENDIF

      PRINTTAB = .FALSE.


C DON'T DO ANYTHING WITH MC...
      IF (COUTAB(FMCKIN) .GT. 0) RETURN

C ONLY USE ON RUNS 22715 - 22718
      IF (ZDSKEY_Nr1 .LT. 22715 .OR. ZDSKEY_Nr1 .GT. 22719) RETURN

      NRCELLS = COUTAB(CALTRU)
      DO I=1,NRCELLS
        CALTRU_ID = I
        CALL GETTAB(CALTRU)
        IF (
     &               CALTRU_CELLNR              .EQ. 11052
     &       .AND.   ABS(CALTRU_IMBAL/CALTRU_E) .GT. 0.8    
     &                                                      )   THEN
            CALL DELTAB(CALTRU)
        ENDIF
      ENDDO

      RETURN
      END
