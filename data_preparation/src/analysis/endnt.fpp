C     ================
      SUBROUTINE ENDNT
C     ================

      IMPLICIT NONE

#include "ezhbook.inc"
#include "zrevt.inc"
#include "common.inc"

      INTEGER icycle

      CALL HCDIR('//ntuple', ' ')
      CALL HROUT(0,icycle, ' ')
      CALL HREND('ntuple')

      print*,' this job ended at run',zrevt_runnr,'event'
     +,zrevt_evtnr(3)

      RETURN
      END

