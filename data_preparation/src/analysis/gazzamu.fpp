 
*     ================
      SUBROUTINE GAZZAMU(rver, igmu)
*     ================
      IMPLICIT NONE
* -------------------------------
*  Author: Gareth Howell - 24/6/97
*
*  INPUT : rver(3)  The primary vertex x,y,z position in centimetres
*
*  OUTPUT: igmu = 0 -> no muon found
*          igmu = 1 -> muon found
*-------------------------------
      REAL Rver(3)
      INTEGER igmu
*
      CALL MULOOP(rver, igmu)
*
      END
