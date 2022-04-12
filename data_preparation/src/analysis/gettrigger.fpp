C
C
C
C     ---------------------------------------------------
      SUBROUTINE GETTRIGGER(IERRT)
C     ---------------------------------------------------
      implicit none
      integer GFLT_var(20)
C     ---------------------------------------------------
C
C  COMMENTS FROM MIKE Wodarczyk
C
C  get_flt_data fills the array GFLT_VAR with the variables necessary 
C    to apply restrictions to GFLT slots 46 and 47 
C    As a side effect, the common block o1cpdt is filled with ALL GFLT 
C    trigger quantities comming from the data (or monte carlo zgana). 
C
C  Output:
C     integer GFLT_VAR(20)  ! array of GFLT quantities needed
C               GFLT_VAR(1) = o1evnt_Subtrg(1) bits 0 to 31   (Needed)
C               GFLT_VAR(2) = o1evnt_Subtrg(2) bite 32 to 63  (Needed)
C               GFLT_VAR(3) = o1subt_Subtrg(1) bits 0 to 31   (Not Needed) 
C               GFLT_VAR(4) = o1subt_Subtrg(2) bite 32 to 63  (Not Needed) 
C               GFLT_VAR(5) = AH_SRTD_s0_ti  srtd timing 0    (Not Needed) 
C               GFLT_VAR(6) = AH_SRTD_s1_ti  srtd timing 1    (Not Needed) 
C               GFLT_VAR(7) = AH_SRTD_s2_ti  srtd timing 2    (Not Needed) 
C               GFLT_VAR(8) = AH_SRTD_s3_ti  srtd timing 3    (Not Needed) 
C               GFLT_VAR(9) = RCAL_EMC_E     CFLT value        (Needed)
C               GFLT_VAR(10) = CAL_E         CFLT E total      (Needed)
C               GFLT_VAR(11) = Rcal_EMCTh    Rcal EMC Threshold bit energy   (Needed)
C               GFLT_VAR(12) = BCAL_EMC_E    Bcal Emc energy   (Needed)
C               GFLT_VAR(13) = CAL_Et        CFLT Et           (Needed)
C               GFLT_VAR(14) = RCAL_ISOE     Sum of Rcal isoe's   (Needed)
C               GFLT_VAR(15:20) = [ not used yet ]
C
C  Implicit output:
C      If you want other GFLT quantities for personal use then
C      place the following includes in your zuanal routine to 
C      access other GFLT variables.
C
C  #include "o1cmn/gflcpd.inc"  
C  C    integer o1cpdt(1000) is the array of ALL GFLT Variables
C  #include "FLT/GFLT/gfltpr.inc"
C  C    parameters in gfltpr allow you to index the array o1cpdt as follows
C  C    REMC = o1cpdt(RCAL_EMC_E)
C
C ======================================================================

#include "common.inc"
c#include "comnt.inc"
#include "fmckin.inc"
#include "partap.inc"
!  /*       <-- this is the necessary ADAMO definitions  */
#include "zrevt.inc"
!  /*        <-- event information  */
#include "zdskey.inc"
#include "gabffw.inc"
#include "vctrhl.inc"
#include "lmresu.inc"
#include "lmposd.inc"
#include "lmeb.inc"
#include "gfltpr.inc"
#include "gfltcr.inc"
#include "ccgsum.inc"
C#include "datatyp.inc"
#include "o1evnt.inc"
!  /*               GFLT Trigger information  */
#include "o1subt.inc"
!  /*               GFLT Subtriggers that fired  */
#include "o2sdec.inc"
#include "o2sbor.inc"
#include "o4sbor.inc"
#include "tltevt.inc"
C#include "gflpar.inc"
c#include "FLT/GFLT/gfltpr.inc"
!  /*               indices for GFLT common block words */
c#include "o1cmn/gflcpd.inc"
#include "gflcpd.inc"
!  /*               for O1CPDT  from o1recon */
c#include "FLT/GFLT/gfltcr.inc"
!  /*               for O1RCEX  */



c Values used by Adi but not be me (yet)
c      integer FLT(4), SLT, TLT(15), DST 
      integer FLT(4), SLT, DST
       integer Ndata,GFLT_Array(999),bcn
c       integer RCAL_EMC_EBP,RCAL_EMC_E,
c     +         RCAL07_ISOL_E
c       PARAMETER ( RCAL07_ISOL_E = 101 )
c       PARAMETER ( RCAL_EMC_E = 96 )
c       PARAMETER ( RCAL_EMC_EBP = 95 )
c       INTEGER RCAL_EBP_STB
c            PARAMETER ( RCAL_EBP_STB = 92 )
c       INTEGER RCAL_E_OUT_BP
c            PARAMETER ( RCAL_E_OUT_BP = 93 )
c      INTEGER RCAL_E_AWAY_B
c            PARAMETER ( RCAL_E_AWAY_B = 94 )
c       INTEGER LUMI_EGAMMA
c            PARAMETER ( LUMI_EGAMMA = 159 )
c       INTEGER LUMI_EE
c            PARAMETER ( LUMI_EE = 160 )


      logical trigmontecarlo
      integer it,ierrt

      LOGICAL FIRST
      data first /.true./

      ierrt=0

      IF (first) then 
       WRITE (*,*) ' ************* '
       WRITE (*,*) '   GETTRIGGER  '
       WRITE (*,*) ' ************* '
       first = .false.
      ENDIF

C>>>>>>>
C
C   Check if data or monte carlo
C
      if (coutab(FmcKin).gt.0) then
         trigmontecarlo = .TRUE.
      else
         trigmontecarlo = .FALSE.
      endif

C Trigger stuff
C
C
C First level trigger
      call vzero(flt,4)
      if(coutab(O1EVNT).gt.0)then
         call FETTAB(O1EVNT,ID,1)     !!! prescaled
         FLT(1) = O1EVNT_Subtrg(1)
         FLT(2) = O1EVNT_Subtrg(2)
      endif
      if(coutab(O1SUBT).gt.0)then
         call FETTAB(O1SUBT,ID,1)     !!! without prescales
         FLT(3) = O1SUBT_Subtrg(1)
         FLT(4) = O1SUBT_Subtrg(2)
      endif
C
C get CFLT data
C
C      call getfltd(imc, cale, bcal_emc, rcal_emc, remcth,
C     &     et, emc_e, tlumi_e, tlumi_g)
C
C      write(6, '(A, I3, 8I8)') 'data',
C     &     imc, cale, bcal_emc, rcal_emc, remcth, et, emc_e,
C     &     tlumi_e, tlumi_g
C
C
C Second level trigger
      SLT = 0
      IF ( COUTAB(O2SDEC).GT.0 ) THEN
         DO IT=1,COUTAB(O2SDEC)
            CALL FETTAB(O2SDEC,ID,IT)
            IF ( O2SDEC_SUBTRIGNO.EQ.2 )
     &           SLT = SLT + IBITS( O2SDEC_SUBTRIGTYPE(1),0,8 )
            IF ( O2SDEC_SUBTRIGNO.EQ.3 ) THEN
               IF ( BTEST(O2SDEC_SUBTRIGTYPE(1),5) ) SLT = SLT + 256
               IF ( BTEST(O2SDEC_SUBTRIGTYPE(1),6) ) SLT = SLT + 512
               IF ( BTEST(O2SDEC_SUBTRIGTYPE(1),7) ) SLT = SLT + 1024
               IF ( BTEST(O2SDEC_SUBTRIGTYPE(1),8) ) SLT = SLT + 2048
            END IF
         END DO
      END IF
C
C
C
C... Third level trigger
      call vzero(TLT,15)
      if(coutab(TLTEVT).gt.0)then
         call FETTAB(TLTEVT,ID,1)
C         TLT(4) = TLTEVT_Subtrg(4)
C         TLT(1) = TLTEVT_Subtrg(1)
         do iT=1,15
          TLT(IT) = TLTEVT_Subtrg(IT)
         enddo
      endif       
C
C... DST bits
      dst = zdskey_tstam11
C
C
C
C.. additional trigger info
C -----------------------------------------------------------------
C  First get the GFLT values for things like RCAL_EMC_E and CAL_E.
C
         If (.not.trigMonteCarlo) then
            CALL O1GTCV(ierrt, O1CPDT)
         else
            call o1vzst(-1)
            call o1ezin(ierrt)
            if (ierrt .ne. 0) then
               write(6, *) 'ERROR getfltd: MC FLT data is missing'
            endif
C     First, just the dead-copy
            MAX_O1CPDT = MAX_O1RCEX
            TOT_O1CPDT = TOT_O1RCEX
            DO  IT = 1, MAX_O1RCEX
               O1CPDT(IT) = O1RCEX(IT)
            ENDDO
            DO  IT = MAX_O1RCEX+1, LIMIT_O1CPDT
               O1CPDT(IT) = -1000
            ENDDO
         endif
         
C
C
C
C
C  Now the common block O1CPDT is filled with GFLT information
C     the indices are given in the file "FLT/GFLT/gfltpr.inc"
C
C
C  There are tons of variables that you might want to save, but 
C    I wrote this specifically to collect the variables used
C    for GFLT(46) in 1996 and 1997.  I will pack them in an 
C    array suitable for passing to the routines that will
C    restrict GFLT(46) to only the parts used for a High Q2
C    F2 analysis.
C
C

         GFLT_VAR(5) = o1cpdt(AH_SRTD_s0_ti)
         GFLT_VAR(6) = o1cpdt(AH_SRTD_s1_ti)
         GFLT_VAR(7) = o1cpdt(AH_SRTD_s2_ti)
         GFLT_VAR(8) = o1cpdt(AH_SRTD_s3_ti)

         GFLT_VAR(9) = o1cpdt(RCAL_EMC_E)
         GFLT_VAR(10) = o1cpdt(CAL_E)
         GFLT_VAR(11) = o1cpdt(RCAL_EBP_STB) + o1cpdt(RCAL_E_OUT_BP)
         GFLT_VAR(12) = o1cpdt(BCAL_EMC_E)
         GFLT_VAR(13) = o1cpdt(CAL_ET)
         GFLT_VAR(14) = o1cpdt(RCAL04_ISOL_E)+o1cpdt(RCAL05_ISOL_E)+
     &        o1cpdt(RCAL06_ISOL_E)+o1cpdt(RCAL07_ISOL_E)

C -----------------------------------------------------------------
C
C
C
C
      call VZERO(TrigDat,15)
      call O1GTCV(Ndata,GFLT_Array)
      TrigDat(1) = Float(GFLT_Array(LUMI_EGAMMA))
      TrigDat(2) = Float(GFLT_Array(LUMI_EE))
C
      Trigdat(3) = Float(GFLT_Array(Rcal_EMC_E))
C
C      TrigDat(2) = Float(GFLT_Array(Rcal_EMC_EBP))
C      TrigDat(3) = Float(GFLT_Array(Rcal_E_OUT_BP))
C      TrigDat(4) = Float(GFLT_Array(Rcal_EBP_STB))
C
      if(coutab(CCGSUM).gt.0) then
         call FETTAB(CCGSUM,ID,1)
         TrigDat(4) = (CCGSUM_ETotalEmc + CCGSUM_ETotalHac)
     +              - (CCGSUM_PzEmc + CCGSUM_PzHac)
      endif
      if(coutab(TLTEVT).gt.0)then
         call FETTAB(TLTEVT,ID,1)
         TrigDat(5) = TLTEVT_Eminpz
      endif
      TrigDat(6) = Float(GFLT_VAR(5))
      TrigDat(7) = Float(GFLT_VAR(6))
      TrigDat(8) = Float(GFLT_VAR(7))
      TrigDat(9) = Float(GFLT_VAR(8))
      TrigDat(10) = Float(GFLT_VAR(9))
      TrigDat(11) = Float(GFLT_VAR(10))
      TrigDat(12) = Float(GFLT_VAR(11))
      TrigDat(13) = Float(GFLT_VAR(12))
      TrigDat(14) = Float(GFLT_VAR(13))
      TrigDat(15) = Float(GFLT_VAR(14))

C>>>>>>>



      RETURN
      END
C


