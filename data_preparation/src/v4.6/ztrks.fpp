      Subroutine z_Trks(Ierr)
      Implicit NONE
*-------------------------------------------------------------------------------
*
*   Called by z_RecGB to fill tracking info (Commons /zRec03/,/zRec04/)
*
*     INPUT:  NONE
*
*     OUTPUT:
*
*       Ierr:   0  <--> Everything is cosher 
*              -1  <--> nTrks in VCTPAR > nT_Max
*
*   Author:   Gennady Briskin, Tel Aviv University
*   Date:     20-Dec-1996
*
*_______________________________________________________________________________
*
#include "partap.inc"
#include "vctrhl.inc"
!  /*    ! Tracking table RAW  */
#include "vcatcal.inc"
!  /*    ! Tracking table RAW at Calorimeter  */
#include "vctpar.inc"
!  /*    ! Tracking table fitted to VERTEX  */
#include "vcparcal.inc"
!  /*    ! Tracking table fitted to VERTEX at Calorimeter  */
*
#include "zdrecgb.inc"
#include "vctvtx.inc"
#include "vcparsec.inc"
#include "vcvtxsec.inc"

*==>  OutPut
      Integer Ierr
      
*==>  Counters
      Integer jTr,I,J
      Logical OK

*==>  Track-Island matching(VCTDCA)
      Real    ptTr,pTrfull
      Real    A1,A3
      Real    terr1,terr2,terr3,terr4

      Real    DistCA,DCATerr,DCACerr
      Real    Coord(3),DCoord(3)

      Logical L_eTrk

c --- additionals for secondary vertex tracks:
      integer  ind_sec, ic1, ic2
      integer  nsvtx, svtxok(nT_max)
      integer  vtx_id
      integer  nT_sec
      integer  hlid(nT_max) 
      real     x_sec, y_sec, z_sec, r_sec

c --- Cut Parameters for BEAM PIPE or CTDINNERWALL
c --- z-vertex range [cm] for secondary vertex tracks
c --- vertices are excluded if (z_sec.lt.ZBPMIN).or.(z_sec.gt.ZBPMAX)
      real       ZBPMIN, ZBPMAX
      parameter (ZBPMIN=-60., ZBPMAX=60.)
c --- approx radial distance of beam pipe [cm]
c --- vertices are excluded if ((r_sec.gt.RBPMIN).and.(r_sec.lt.RBPMAX))
      real       RBPMIN, RBPMAX
      parameter (RBPMIN=5., RBPMAX=9.)
c --- approx radial distance of CTD inner wall [cm] 
c --- vertices are excluded if ((r_sec.gt.CTDIWMIN).and.(r_sec.lt.CTDIWMAX))
      real       CTDIWMIN, CTDIWMAX
      parameter (CTDIWMIN=14.,CTDIWMAX=18.)

      logical  secfirst, rnotok
      data secfirst / .true. /

*==>  Ini
      Ierr    = 0

      eTrMat  = 0
      eTrNear = 0
      eTrDCA  = 999.0

*==>  Do we have enough space for tracks ??
      nT      = CouTab(VCTPAR)
      If (nT.GT.nT_Max) Then
       Write(6,*) ' '
       Write(6,*) '*==> z_Trks: more than',nT_Max,' tracks ??'
       Write(6,*) ' '
       Ierr = -1
       Return
      EndIf

*==>  Do we have tracks in this event ??
      If (nT.EQ.0) Return

*==>  get eDIS Track Info.........................
      If (eNdis.GT.0) Then
       
       Coord( 1) = eXdis
       Coord( 2) = eYdis
       Coord( 3) = eZdis
       DCoord(1) = 5.0
       DCoord(2) = 5.0
       DCoord(3) = 5.0
        
*==>  Do the match :~)
       Do I=1,nT
        Call VCTDCA(I,Coord,DCoord,DistCA,DCATerr,DCACerr)
        If (DistCA.GE.0..AND.DistCA.LE.eCut_DCA) Then
                eTrMat  = eTrMat + 1
         vcteID(eTrMat) = I
         If (DistCA.LT.eTrDCA) Then
          eTrDCA  = DistCA
          eTrNear = I
         EndIf
        EndIf
       EndDo
       
      EndIf
*
*==>  Hadronic Track Info
      jTr = 0
      Do 100 I=1,nT
*
       Do J=1,eTrMat
        If(I.EQ.vcteID(J)) GOTO 100
       EndDo
*     
       Call FetTab(VCTPAR  ,ID,I)
       Call FetTab(VCPARCAL,ID,I)
       Call NatRel(VCTPAR,VCTPAR_VCTRHL,VCTRHL,OK)
       If ( .NOT. OK ) Then
        Print *,'z_Trks: No such VCTRHL track ',I
        STOP
       End If         
*     
*==>>  Calculate momentum of the track in the Lab (ZEUS-Note 96-013)
*     
       ptTr    = VCTRHL_pgevc/Sqrt(1.0+VCTRHL_tdip**2)
       ptTr    = ptTr*ABS(VCTRHL_qovr/VCTPAR_Par(3))
       pTrfull = ptTr/Sin(VCTPAR_Par(1))
       
*==>   Calculate momentum error on the track
       A1      = VCTPAR_Par(1) 
       A3      = VCTPAR_Par(3)
       terr1   = 1.0/TAN(A1)
       terr2   = 1.0/A3
       terr3   = 2.0*terr1*terr2
       terr4   = pTrfull*Sqrt(VCTPAR_Cov(1)*terr1**2 +
     +                        VCTPAR_Cov(6)*terr2**2 +
     +                        VCTPAR_Cov(4)*terr3    
     +                       )
     
*==>   Track variable 
       jTr         = jTr + 1
       vcthID(jTr) = I
       hlid(  jTr) = vctpar_VCTRHL
       pTr(   jTr) = pTrfull 
       thTr(  jTr) = VCTPAR_Par(1) 
       phTr(  jTr) = VCTPAR_Par(2)
       dpTr(  jTr) = terr4
       swmTr( jTr) = VCPARCAL_KODSWM
*
C ---  AQ modification (trust tracking more in transition region) ---
       IF(mod(swmTr(jTr)/100,10).ge.8) Then
         ptMax_CUT = ptMax_CUT9
       ELSE
         ptMax_CUT = ptMax_CUT7
       ENDIF

       If(mod(swmTr(jTr)/100,10).GE.sl_CUT    .AND.
     &    ptTr                  .GE.ptMin_CUT .AND.
     &    ptTr                  .LE.ptMax_CUT )Then
        qltTr(jTr) = 1
       ELSEIf(mod(swmTr(jTr)/100,10).GE.sl_CUT     .AND.
     &        ptTr                  .GE.ptMin_CUT  .AND.
     &        ptTr                  .LE.muon_Ptmax .AND.
     &        AQ_options                               )Then
        qltTr(jTr) = 2
       Else
        qltTr(jTr) = 0
       EndIf

*
 100  Continue
      nT = jTr
*
c     ===============================================================
c --- if you want also secondary vertex tracks: 
      if (sec_vtx_yes) then
         if (secfirst) then
            secfirst = .false.
            ind_sec = getind(VCPARSEC,'VCTRHL')
         endif
c ---  Return if there aren't any secondary vertices
c         if (Coutab(VCPARSEC).EQ.0) RETURN
         if (Coutab(VCPARSEC).EQ.0) goto 222
c --- clear vertex array
         nsvtx = 0
         do i=1,nT_max
            svtxok(i) = 0
         enddo
            
c --- look for electron track here if hasn't been found in VCTPAR:
c         If ((eNdis.GT.0).and.(eTrNear.eq.0)) Then

c --- which of the secondary vertices are accepted? 
         do i=1,Coutab(VCVTXSEC)
            call fettab(VCVTXSEC,id,i)            
            x_sec = vcvtxsec_v(1)
            y_sec = vcvtxsec_v(2)
            z_sec = vcvtxsec_v(3)
            r_sec = sqrt(x_sec*x_sec + y_sec*y_sec)
            rnotok = .true.
            if ((r_sec.gt.RBPMIN).and.(r_sec.lt.RBPMAX)) then
               rnotok = .false.
            endif
            if ((r_sec.gt.CTDIWMIN).and.(r_sec.lt.CTDIWMAX)) then
               rnotok = .false.
            endif
            if ((z_sec.lt.ZBPMIN).or.(z_sec.gt.ZBPMAX).or.(rnotok)) then
               nsvtx = nsvtx + 1
               svtxok(nsvtx) = VCVTXSEC_id
            endif
         enddo

c --- grab the tracks (go back to VCTRHL)
         nT_sec = Coutab(VCTRHL)
         Do I=1,nT_sec
            call fettab(VCTRHL,id,i)
            call fettab(VCATCAL,id,i)
            call nafrel(VCTRHL,VCPARSEC_VCTRHL,VCPARSEC,ind_sec,ic1,ic2)
            if (ic1.eq.ic2) then
               call fettab(VCPARSEC,ind_sec,ic1)
               vtx_id = VCPARSEC_ProducedAt
c --- check if vertex not ok, if not goto 200 , ifyes continue               
               do j=1,nsvtx
                  if (vtx_id.eq.svtxok(j)) goto 200
               enddo
c --- check if VCTRHL-track no not in VCTPAR, if so, do not take 
c        this track but print out error message
               do j=1,nT
                 if (i.eq.hlid(j)) then
                     print*, 'z_Trks> MESSAGE'
                     print*, 'same VCTRHL track pointing to primary ',
     & 'and secondary vertex' 
                     goto 200
                  endif
               enddo
*     
**-->> Calculate momentum of the track in the Lab (ZEUS-Note 96-013)
*     
               ptTr    = VCTRHL_pgevc/Sqrt(1.0+VCTRHL_tdip**2)
               ptTr    = ptTr*ABS(VCTRHL_qovr/VCPARSEC_Par(3))
               pTrfull = ptTr/Sin(VCPARSEC_Par(1))
               
*-->   Calculate momentum error on the track
               A1      = VCPARSEC_Par(1) 
               A3      = VCPARSEC_Par(3)
               terr1   = 1.0/TAN(A1)
               terr2   = 1.0/A3
               terr3   = 2.0*terr1*terr2
               terr4   = pTrfull*Sqrt(VCPARSEC_Cov(1)*terr1**2 +
     +              VCPARSEC_Cov(6)*terr2**2 +
     +              VCPARSEC_Cov(4)*terr3    
     +              )
     
*-->   Track variable 
               jTr         = jTr + 1
               if (jTr.gt.nT_max) goto 222
               vcthID(jTr) =  - i
               pTr(   jTr) = pTrfull 
               thTr(  jTr) = VCPARSEC_Par(1) 
               phTr(  jTr) = VCPARSEC_Par(2)
               dpTr(  jTr) = terr4
               swmTr( jTr) = VCATCAL_KODSWM
*
               If(mod(swmTr(jTr)/100,10).GE.sl_CUT.AND.
     &              ptTr.GE.ptMin_CUT.AND.ptTr.LE.ptMax_CUT )Then
                  qltTr(jTr) = 1
               Else
                  qltTr(jTr) = 0
               EndIf
            endif
*     
 200        continue
         enddo
      endif
 222  continue

      if (printflag) then
         print*, 'VCTRHL table:'
         write(*,*) 'ID ,Q,SWM,  P      ,  THE    ,  X,Y,Z on CAL -->'
         Do I=1,Coutab(VCTRHL)
            call fettab(VCTRHL,id,i)
            call fettab(VCATCAL,id,i)
            write(*,231) VCTRHL_ID, vcatcal_q, vcatcal_KODSWM, 
     &           vctrhl_pgevc, atan(vctrhl_tdip), 
     &           vcatcal_x, vcatcal_y, vcatcal_z
         enddo
         
         print*, ' '
         print*, 'Track SELECTION'
         write(*,*) 'ID ,Q,SWM,  P      ,  THE    ,  PHI '
         do i=1,nT
         write(*,230) hlid(i),qltTr(i),swmTr(i),pTr(i),thTr(i),phTr(i)
         enddo
         do i=nT+1,jTr
         write(*,230) vcthID(i),qltTr(i),swmTr(i),pTr(i),thTr(i),phTr(i)
         enddo
      endif
      
      if (jTr.gt.nT_max) then
         print*, 'z_Trks> too many tracks in event: '
         print*, 'z_Trks> only primary vertex Tracks taken '
      else
         nT = jTr
      endif
 230  format (i3,1x,i1,1x,i3,3(f9.4,1x))
 231  format (i3,1x,i2,1x,i3,5(f9.4,1x))
      End

