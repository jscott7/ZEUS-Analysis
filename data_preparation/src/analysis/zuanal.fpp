C     =========================
      SUBROUTINE ZUANAL (IWANT)
C     =========================
C
C     ----------------------------------------------------------------
C     Routine Called Every Event from EAZE.
C     Part of the Initilisation is done here and not in ZUINIT.
C --------------------------------------------------------------------

      IMPLICIT NONE

#include "partap.inc"
#include "zrevt.inc"
#include "o1setu.inc"
#include "fmckin.inc"
#include "zdskey.inc"
c#include "ctrlcm.inc"
#include "tempcm.inc"
#include "ezhbook.inc"
#include "ccgsum.inc"
#include "common.inc"
#include "caltru.inc"
#include "cconsa.inc"

      INTEGER      IERR,IWANT
      LOGICAL      FIRST
      DATA         FIRST / .TRUE./
      real vtx(3)
      integer i, j
      logical MC

c For PCCnds
      real energycuts(10)
c      data energycuts /0.08,0.150,-0.1,-0.2,-0.2,-2.0,-4.0,-8.0,
c     &     -9.0,-10.0/
      data energycuts /-1.0,-1.0,-1.0,-1.0,-1.0,
     &                 -1.0,-1.0,-1.0,-1.0,-1.0/
      integer nrcond, ierrcond

c For noise96    
      real ethresh(2),imbacut
      integer noisycells, ierrnoise
      logical islfl,delfl
      data ethresh /0.10,0.16/
      data imbacut /0.7/
      integer nsigma
      data nsigma / 3 /
      data islfl /.false./
      data delfl /.true./

c VCEAZE stuff 
      integer jinit, jtrack, jfit, jdedx, jvertex, jerr
      logical Lfirst
      data Lfirst  /.true./  

      IWANT    = 0              !!! Do not Write-Out Event
      do i = 1, 3
         vtx(i) = 0.
      end do
C --------------------------
C --- Init on first call ---
C --------------------------
      IF (FIRST) THEN

         WRITE(*,*)
         WRITE(*,*) ' Basic initialisation'
         WRITE(*,*)

c ------------------------------------------------
c CCGEOM, SRGEOM, get year, Check if real data/MC.
c ------------------------------------------------         
         CALL BAS_INIT(MC)
         FIRST = .FALSE.

         write(*,*) 'Calling VCEAZE for multiple vertexing'
         call VCEAZE(1,0,0,0,1,ierr)
         If (ierr.ne.0) write(*,*) 'tering-VCEAZE --- ',ierr 

         Call SREAZE(-1,0,1,ierr)
         If (ierr.ne.0) write(*,*) 'tering-SREAZE --- (init)',ierr

      ENDIF

      montecarlo = MC
c --------------------------------------
c Now for event by event processing
c --------------------------------------
c Reject events not selected by EVTAKE 
c --------------------------------------
      if (.not.montecarlo) then
         call evtake(ierr)
         if (ierr.eq.0) then
            write(*,*),'Event rejected by EVTAKE'
            return
         end if
      end if

C --------------------------
C --- Run & event number ---
C --------------------------
      if (coutab(ZREVT).ne.0) then
         CALL FETTAB (ZREVT,ID,1)
         RUN_NUM = ZREVT_RunNr
         EV_NO = ZREVT_EvtNr(3)
         FLT_NO = ZREVT_EvtNr(1)
      else
         write(*,*) 'ZUANAL : ZREVT empty - Abort event'
         return
      end if

c --------------------------------------
c SRTD reconstruction
c --------------------------------------
      call sreaze(0,0,1,ierr)
      if (ierr.ne.0) write(*,*) 'ERROR IN SREAZE'

      CALL CELSWP               !! Correct miscabled BCAL EMC Cells
      CALL TRKFIL               !! Tracking stuff
      CALL SPRKREJ		!! Spark Rejection
      call rmspark(Sparks)

c ---------------------------------------------
c Run VCEAZE for secondary vertexing for use by
c ZufosNT
c ---------------------------------------------
c      call preptracks(ierr)

c -----------------------
c Get Trigger Information
c -----------------------
c FLT
c -----------------------
      call gettrigger(ierr)


      SLT_EMPZ = 0.

c ----------------------------------------------
c Get Beam Energies from O4SBOR table or FMCKIN
c ----------------------------------------------
      CALL BEAMFIL             

      CALL GETLUMI               !! Lumi Detector Info

c      do i = 1, coutab(caltru)
c         Caltru_ID = i
c         call gettab(caltru)
c         caltru_cconsa=INULL
c         call reptab(caltru)
c      end do

c --------------------------------------------------------  
c Do Calorimeter energy corrections
c --------------------------------------------------------
c      call rcalcorr(ierr)       !! Old Routine
      call rcalcorr6(ierr)      !! RcalCor 6 (Recommended)

c--------------------------------------------------------
c Do Noise Routine
c --------------------------------------------------------     
      call pccnds(energycuts,nrcond,ierrcond)
      if (ierrcond.ne.0) write(*,*)'PCCnds:error'  
      call NOISE96S(run_num,ethresh,imbacut,nsigma,
     &     islfl,delfl,noisycells,ierrnoise)
      if (ierrnoise.ne.0) print*,'error in noise96', ierr
      call pccnds(energycuts,nrcond,ierrcond)

c --------------------------------------------------------
c Run Sinistra95 & FINDIS option 1
c --------------------------------------------------------
      call setvtx(vtx)

      CALL GETSINISTRA(vtx, year)           

c Reject events where no electron is found.
      if (CAL_EE.lt.0.) then
         write(*,*) '**** NO SINISTRA ELECTRON FOUND! ****'
         return
      end if
      call setvtx(vtx)
c ------------------------
c Run ZufosNT
c ------------------------
      CALL GETZUFOSNT(vtx)

c ------------------------
c Get Presampler Info.
c ------------------------
      CALL GETPRES(vtx)

      CALL GETCAL(vtx)               !! Calorimetry stuff

      if (montecarlo) then		
         CALL MR2FIL            !! Fill MC true
      end if
    
      call bgdfil(montecarlo)
      CALL HFNT(10)             !! Fill ntuple

      RETURN
      END
      






