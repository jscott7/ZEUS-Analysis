C       ===============================================
        SUBROUTINE FINDIS95( Iopt, OutCut, Cand, Ierr )
C       ===============================================
C        
C       Author       Ralph Sinkus
C       Modified     G. Briskin 
c                    Option 4 added
c       Modified     A. Lopez-Duran, O.Ruske, S. Schlensted
c                    Found bug in option 4, and fixed it
c       Modified     A.Lopez-Duran
c                    Option 5 added
c
C...    Routine for siNISTra
C...    Finds among candidates stored in CANDAT the DIS e-
c
C       INPUT :
c
c...    Iopt = 1  <-> candidate with highest probability from siNISTra.
c...                  If 2 or more candidates have the same probability,
c...                  the one with highest energy is going to be picked up
c...                  as the DIS e-
c...                  !!! Candidates with Yel>0.99 are NOT considered. !!!!
c
c...         = 2  <-> candidate with highest E
c
c...         = 3  <-> candidate with highest Pt
c
c...         = 4  <-> candidate with highest probability from siNISTra.
c...                  If 2 or more candidates have the same probability,
c...                  the one with highest energy is going to be picked up
c...                  as the DIS e-. (doVCTDCA has to be set .TRUE. in SIDAT95)
c...                  !!! Candidates with Yel>0.99 are NOT considered. !!!!
c...                  !!! Candidates in BCAL, FCAL(Re>70) and RCAL(Re>50) which
c...                      DO NOT have a track pointing to it are NOT considered !!!
c...
c...         = 5  <-> same as option 4, but track region calculated using the 
c...                  event vertex instead of z vertex = 0.0   
c
c       OutCut = 0.9   <> Cut in probability of NN. This value is recommended
c
c
c       OUTPUT :
c
c...   Cand           <> Number of candidate which is supposed to be the DIS e-
c...    Ierr = 0       <> No error
c...         = 1       <> WARNING !!!!!!!!!!!!!!!!
c...
#include "partap.inc"
#include "sidat95.inc"
#include "vctvtx.inc"

      INTEGER Iopt,Cand,Ierr
      REAL    OutMax, Emax, PtMax, Pt, OutCut
      REAL    Yel
      INTEGER I,K,NoutMax
      logical first 
      data    first /.true./
*-->  initial energy of e+ and other stuff
      Real       eE0
      Parameter( eE0=27.52 )

      Real       twoeE0
      PARAMETER( twoeE0 = 2.0*eE0 )

      Real       eE,eY,eX,eZ,eR
      Real       Zvtx
      Logical    LGOOD
*
* Variables for option 5
* -------------------------------------------------
      real       vz, d_FCal, d_RCal, r_FCal, r_RCal
      real       line, point
*
*                 3 super layers  4 super layers
      real        t_slay,         f_slay,
*                 f side of CTD   r side of CTD
     &            f_ctd ,         r_ctd ,
*                 distance to FCal, RCal 
     &            f_cal ,         r_cal ,            
*                 radius of BCal
     &            r_BCal,
     &            f_rat ,         r_rat
*
      parameter ( t_slay = 35.,   f_slay = 45.4, 
     &            f_ctd  = 106.,  r_ctd  = 101.,
     &            f_cal  = 226.1, r_cal  = 152.1,
     &            r_BCal = 125.6)
* -------------------------------------------------
*
c...    Initialize
      Cand    = 0
      Ierr    = 0
      NoutMin = 0
      Emax    = 0.0
      PtMax   = 0.0

      if(first)then
         first =.false.
       WRITE(6,*) '**********************************'
       WRITE(6,*) '* This is findis95   Version 2.0 *'
       WRITE(6,*) '*   This version includes        *'
       WRITE(6,*) '*   options 4 and 5              *'
       WRITE(6,*) '**********************************'
      endif
C...    check for negative OutCut
      IF ( OutCut.LT.0.0 ) THEN
       WRITE(6,*) '**********************************'
       WRITE(6,*) '* This is siNISTra95 Version 1.1 *'
       WRITE(6,*) '*   Probability definition has   *'
       WRITE(6,*) '*   changed !                    *'
       WRITE(6,*) '*   Electromagnetic objects      *'
       WRITE(6,*) '*   populate now the probability *'
       WRITE(6,*) '*   region close to one !        *'
       WRITE(6,*) '*   Chose OutCut = 0.9 !         *'
       WRITE(6,*) '* Program stopped by siNISTra95  *'
       WRITE(6,*) '**********************************'
       STOP '---- Program aborted by siNISTra95 ----'
      EndIf



C...    No candidates available
      If ( Ncand.EQ.0 ) GOTO 1000

*-->  VC - Reconstructed Vertex....................
      If ( CouTab(VCTVTX).GT.0 ) Then
       Call FetTab(VCTVTX,ID,1)
       Zvtx = VCTVTX_V(3)
      Else
       Zvtx = 0.0
      EndIf

*-->  Check if the Vertex is cosher, a la elecpo
      If (Zvtx.LT.-148. .OR. Zvtx.GT.235.) Then
       Zvtx = 0.0
      EndIf

      IF ( Iopt.EQ.1 ) THEN
C...     find highest probability among candidates with OUT values above OutCut
       OutMax = -1.0
       Do I=1,Ncand
        eE  = CANDAT(2,I)
        eX  = CANDAT(3,I)
        eY  = CANDAT(4,I)
        eZ  = CANDAT(5,I) - Zvtx
        Yel = 1. -  (eE/twoeE0) *
     &       (1. - eZ/Sqrt(eX**2+eY**2+eZ**2))

        IF ( CANDAT(1,I) .GE. OutCut .AND.
     &       Yel.LT.0.99             .AND.
     &       CANDAT(1,I) .GT. OutMax
     &       ) OutMax = CANDAT(1,I)
       EndDo
C...      found at least a candidate ?
       IF ( OutMax.EQ.-1.0 ) RETURN


C...      how often does this maximal out values appear ?
       DO I=1,Ncand
        IF ( CANDAT(1,I) .EQ. Outmax ) NoutMax = NoutMax + 1
       END DO
       IF ( NoutMax .EQ. 1 ) THEN
        DO I=1,Ncand
         IF ( CANDAT(1,I) .EQ. OutMax ) Cand = I
        END DO
       ELSE
C...        find among those candidates the one with highest E
        DO I=1,Ncand
         IF ( OutMax .EQ. CANDAT(1,I) .AND.
     &        Emax   .LT. CANDAT(2,I)
     &        ) Emax = CANDAT(2,I)
        END DO
        DO I=1,Ncand
         IF ( CANDAT(2,I) .EQ. Emax ) Cand = I
        END DO
       END IF


      ELSE IF ( Iopt .EQ. 2 ) THEN
       DO I=1,Ncand
        IF ( CANDAT(1,I).GE.OutCut .AND.
     &       Emax .LT. CANDAT(2,I)
     &       ) THEN
         Emax = CANDAT(2,I)
         Cand = I
        END IF
       END DO


      ElseIf ( Iopt .EQ. 3 ) Then

       DO I=1,Ncand
        IF ( CANDAT(1,I).GE.OutCut ) THEN
         Pt = CANDAT(2,I) * ( (CANDAT(3,I)**2+CANDAT(4,I)**2) /
     &        (CANDAT(3,I)**2+CANDAT(4,I)**2+CANDAT(5,I)**2)
     &        )**0.5
         IF ( Ptmax .LT. Pt ) THEN
          Ptmax = Pt
          Cand  = I
         END IF
        END IF
       END DO

      ElseIf (Iopt.EQ.4.AND.doVCTDCA) Then

C...      find highest probability among candidates with OUT values above OutCut
       OutMax = -1.0
       Do I=1,Ncand

        eE  = CANDAT(2,I)
        eX  = CANDAT(3,I)
        eY  = CANDAT(4,I)
        eZ  = CANDAT(5,I)
        eR  = Sqrt(eX**2+eY**2)

        LGOOD = .TRUE.

        If ( CANDAT(15,I).EQ.0 ) Then
           If( eZ.LT.0. ) then
              if( eR.GT.50. ) LGOOD=.FALSE.
           else
              if( eR.GT.70. ) LGOOD=.FALSE.
           endif
        endif

        eZ  = eZ - Zvtx
        Yel = 1. -  (eE/twoeE0) *
     &       (1. - eZ/Sqrt(eX**2+eY**2+eZ**2))

        If (Yel.LT.0.99.AND.LGOOD) Then
         If(CANDAT(1,I).GE.OutCut.AND.CANDAT(1,I).GT.OutMax)Then
          OutMax = CANDAT(1,I)
         EndIf
        EndIf

       EndDo

C...      found at least a candidate ?
       If ( OutMax.EQ.-1.0 ) Return

C...      how often does this maximal out values appear ?
       Do I=1,Ncand
        If (CANDAT(1,I).EQ.Outmax) NoutMax = NoutMax + 1
       EndDo

       If ( NoutMax .EQ. 1 ) Then

        Do I=1,Ncand
         If (CANDAT(1,I).EQ.OutMax) Then
          Cand = I
          Return
         EndIf
        EndDo

       Else
*==>    Find among those candidates the one with highest E
        Do I=1,Ncand
         If(CANDAT(1,I).EQ.OutMax.AND.Emax.LT.CANDAT(2,I))Then
          Emax = CANDAT(2,I)
          Cand = I
         EndIf
        EndDo

       EndIf

         
      ElseIf ( Iopt.EQ.5.AND.doVCTDCA) Then
         
C... FINDIS option 5
C
         OutMax = -1.0
         Do I=1,Ncand
* Cand position            
            eE  = CANDAT(2,I)
            eX  = CANDAT(3,I)
            eY  = CANDAT(4,I)
            eZ  = CANDAT(5,I)
            eR  = Sqrt(eX**2+eY**2)
            vz  = zvtx
* Correct distance to cal. with vertex            
            d_FCal = f_cal - vz
            d_Rcal = r_cal + vz
* Ini. flag            
            LGOOD  = .TRUE.
            If ( CANDAT(15,I).EQ.0 ) Then
* If cand. no track
               If     ( eZ.LT.-140. ) then                 !e in RCal
* Calculate track region for RCal
                  r_rat  = f_slay / (r_ctd + vz)
                  r_Rcal = d_Rcal * r_rat
* Look if vtx inside CTD and electron in track region
                  if( eR.GT. r_RCal .and.
     &                vz.gt.-r_ctd) LGOOD=.FALSE.
               elseif ( eZ.GT. 220. ) then                 !e in FCal
* Calculate track region for FCal
                  f_rat  = f_slay / (f_ctd - vz)
                  r_Fcal = d_Fcal * f_rat
* Look if vtx inside CTD and electron in track region
                  if( eR.GT. r_FCal .and.
     &                 vz.lt. f_ctd) LGOOD=.FALSE.
               else                                        !e in BCal
                  if( vz.lt.f_ctd.and.vz.gt.-r_ctd)then    !vtx inside CTD
* If electron in BCal and electron position and vtx in CTD region, always track
                     if( eZ.lt.f_ctd.and.eZ.gt.-r_ctd)then !z inside CTD
                        LGOOD = .FALSE.
                     else                                  !z outside CTD
                        if(eZ.lt.vz)then                   !z position wrt vtx
* Calculate intercept point between vtx-electron and CTD
                           line   = R_ctd
                           point  = r_BCal * ( line + vz) / (vz - eZ) 
* If intercept point greater than 5 superlayers, electron in CTD region
                           if (point.gt.f_slay) LGOOD = .FALSE.
                        else
                           line   = f_ctd
                           point  = r_BCal * ( line - vz) / (eZ - vZ) 
                           if (point.gt.f_slay) LGOOD = .FALSE.
                        endif
                     endif
                  else                                      !vtx outside CTD
                     write(6,*)'This is findis95 option 5'
                     write(6,*)'WARNING!!!, vtx outside CTD'
                     write(6,*)'Do not trust the result!!!!'
                     if(vz.gt.f_ctd.and.vz.gt.eZ)then
                        line   = f_ctd
                        point  = r_BCal * ( line + vz) / (vz - eZ) 
                        if (point.lt.f_slay) LGOOD = .FALSE.
                     elseif(vz.lt.-r_ctd.and.vz.lt.eZ)then
                        line   = r_ctd
                        point  = r_BCal * ( line - vz) / (eZ - vz) 
                        if (point.lt.f_slay) LGOOD = .FALSE.
                     endif
                  endif
               endif
            endif
* Continue as other options
        eZ  = eZ - Zvtx
        Yel = 1. -  (eE/twoeE0) *
     &       (1. - eZ/Sqrt(eX**2+eY**2+eZ**2))

        If (Yel.LT.0.99.AND.LGOOD) Then
         If(CANDAT(1,I).GE.OutCut.AND.CANDAT(1,I).GT.OutMax)Then
          OutMax = CANDAT(1,I)
         EndIf
        EndIf

       EndDo

C...      found at least a candidate ?
       If ( OutMax.EQ.-1.0 ) Return

C...      how often does this maximal out values appear ?
       Do I=1,Ncand
        If (CANDAT(1,I).EQ.Outmax) NoutMax = NoutMax + 1
       EndDo

       If ( NoutMax .EQ. 1 ) Then

        Do I=1,Ncand
         If (CANDAT(1,I).EQ.OutMax) Then
          Cand = I
          Return
         EndIf
        EndDo

       Else
*==>    Find among those candidates the one with highest E
        Do I=1,Ncand
         If(CANDAT(1,I).EQ.OutMax.AND.Emax.LT.CANDAT(2,I))Then
          Emax = CANDAT(2,I)
          Cand = I
         EndIf
        EndDo

       EndIf
         
         
      Else
       Write(6,*) 'FINDIS95: Option not yet implemented'
      EndIf

 1000 CONTINUE

      RETURN
      END


