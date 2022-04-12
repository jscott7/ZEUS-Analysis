      FUNCTION ISTUD(KID,DISMIN,GENKEY
     +              ,ICHA,IPAPRT,IPACOD,IDELSA )
C# ----------------------------------------------------------------
C#
C@#   ISTUD : select stable final state particles from FMCKIN (MC TRUTH)
C
C#     FUNCTION ISTUD(KID,DISMIN,GENKEY,ICHA,IPAPRT,IPACOD,IDELSA )
C#    ==============================================================
C#
C#     PURPOSE :  scans FMCKIN and select stable final state particles
C#                a particle is regarded as a "stable" particle if :
C#          a) on ZDIS level :
C#               it has ISTHEP flag = 1
C#          b) on MOZART/ZEPHYR level :
C#               it has been processed by GEANT
C#              BUT NOT if :
C#                it is created in some secondary interaction in the UCAL
C#                     or somewhere outside the CTD
C#            In case of decay in flight ( simulated by GEANT)
C#             the daughter particles are kept if the decay is close
C#             to the I.P. ( < DISMIN )
C#             else the mother particle is kept and the daughters are
C#              not counted
C#
C#    ENVIRONMENT required :
C#
C#      a)  needs ADAMO tables
C#      b)  FMCKIN table MUST be scanned INVERSE ! i.e :
C#
C#       DO <label> KID=COUTAB(FMCKIN),1,-1
C#         icheck= ISTUD(KID ....)
C#
C#
C#    INPUT   :  KID = row no. in FMCKIN (FMCKIN_ID)
C#             DISMIN=  parameter , min. distance of decay vertex required
C#                     if distance < DISMIN  daughter particle declared
C#                     as stable, mother is not taken
C#
C#
C#   OUTPUT  :  ISTUD = 1 final state particle seen in the apparatus
C#                    = 0 : else
C#              ICHA  = -1/0/1 for pos./neutral/neg. charge
C#              IPAPRT=  no. line in FMCPRT
C#              IPACOD=  LUND PART. CODE
C#              IDELSA=  FMCKIN_ID of scattered electron
C#                    =  0 :  else
C#
C#   error return :  ( "ill" tracks )
C#         ISTUD = -777 : no mother for decay paricle found
C#         ISTUD = -888 : particle has no entry in FMCVTX
C#         ISTUD = -999 : particle has no entry in FMCPRT
C#
C#    AUTHOR  :  Nikolaj Pavel                                  29.2.92
C#
C#    modification log :
C#     13.7.92  : bug fix : IDELSA # 0 only for scatt. electron also for
C#                          HCL,HAR ( i.e. MC with rad. photons)
C#    21.1.94  : NPA    : for newest HL6 ( heracles+lepto6.1)
C#                        FMCKIN_daughter for scatt. elec wrong filled
C#                        ==> IDELSA = 0 always returned
C#                        ===> way out : take IDELSC from /LEPKIN/ which is
C#                        filled by VGKIN
C#
C# ----------------------------------------------------------------
*
      IMPLICIT NONE
*
#include "partap.inc"
#include "fmckin.inc"
#include "fmcvtx.inc"
#include "fmcpcd.inc"
#include "fmcprt.inc"
*
      CHARACTER*4 GENKEY
      INTEGER ISTUD,ICHA,IPACOD,IDELSA,IPAPRT
      INTEGER ISTGEN ,KID
      REAL DISMIN ,VMOD
      LOGICAL OK
*
      ISTUD = 0
      ICHA = 0
      IPACOD = 0
      IDELSA = 0
*
      CALL FETTAB(FMCKIN,ID,KID)
C
      IF( FMCKIN_ISTHEP .LT.1000) THEN
C ----  NOT  processed by GEANT or on ZDIS level
C
        IF(FMCKIN_ISTHEP.EQ.1 .AND. .NOT.FMCKIN_DECAY) THEN
C
C ---  it's a stable one -> get particle code/charge
          GOTO 1000
        ELSE
C
C ---       not stable
          GOTO 8900
        ENDIF
C
      ELSE
C
        ISTGEN = MOD(FMCKIN_ISTHEP,1000)
C
        IF(ISTGEN.EQ.0) THEN
C
C ---   produced by GEANT
C
          IF(FMCKIN_DECAY) THEN
C ---     not stable
            GOTO 8900
          ELSE
            CALL NATREL(FMCKIN,FMCKin_PRoducedAt,FMCVTX,OK)
            IF(OK) THEN
              IF(FMCVTX_TYPE .NE. 'DCAY') THEN
C ---    reject particles produced in other process than decay processes
C         in order to avoid double counting
                GOTO 8900
C
              ELSE
C ---    cut on distance from nominal I.P.
                IF(VMOD(FMCVTX_R,3) .LT. DISMIN ) THEN
C ---         stable decay products from decay close to I.P.
C                take it
                  GOTO 1000
                ELSE
C ---         decay to far away => regard mother as stable
                  CALL NATREL(FMCKIN,FMCKin_DaughterOf,FMCKIN,OK)
                  IF(OK) THEN
                    FMCKIN_DECAY = .FALSE.
                    CALL REPTAB(FMCKIN)
                    FMCKIN_ID = KID
                    CALL GETTAB(FMCKIN)
                  ELSE
                    WRITE(6,'(A,I4)')
     +     ' **** WARNING NO mother for decay particle FMCKIN_ID='
     +              ,FMCKIN_ID
                    ISTUD = -777
                  ENDIF
                  GOTO 8900
                ENDIF
              ENDIF
            ELSE
C ---      "ill" track
              WRITE(6,'(A/A,I4,A,I6)')
     +   '**** WARNING FROM >ISTUD< : track without production vertex'
     +  ,'      *** FMCKIN_ID = ',ID,'  FMCKIN_FMCPRT = ',FMCKIN_FMCPRT
              ISTUD = -888
              GOTO 8900
            ENDIF
          ENDIF
        ELSE
C ---    particle from ZDIS touched by GEANT
        IF(ISTGEN.EQ.1 .AND. .NOT.FMCKIN_DECAY) THEN
C
C ---  it's a stable one -> get particle code/charge
            GOTO 1000
          ELSE
C
C ---       not stable
            GOTO 8900
          ENDIF
        ENDIF
C
      ENDIF
C
 1000 CONTINUE
C
C --- particle stable , so take it
      ISTUD = 1
      CALL NATREL(FMCKIN,FMCKin_FMCPrt,FMCPRT,OK)
      IF(OK) THEN
        IPAPRT = FMCPRT_ID
        ICHA = FMCPRT_CHARGE
        FMCPCD_ID = FMCPRT_ID
        CALL GETTAB(FMCPCD)
        IPACOD = FMCPCD_LUNDCOD
        IF(    FMCKIN_DAUGHTEROF.EQ.1
     +   .and. (ABS(IPACOD).EQ.7.OR. ABS(IPACOD).EQ.8) )
     +         IDELSA = FMCKIN_ID
        RETURN
      ELSE
        WRITE(6,'(A,I5)')
     + ' **** WARNING no entry in FMCPRT table for particle FMCKin_ID='
     + ,FMCKIN_ID
        ISTUD = -999
        RETURN
      ENDIF
C
 8900 CONTINUE
C --- reject unstable
C
      RETURN
      END
