C     =============================================
      SUBROUTINE CCCXYZ (CellNr, Cx, Cy, Cz, Error)
C     =============================================
      IMPLICIT NONE
C
C ************************************************************************
C
C PURPOSE:       Calculation of cell centers from offline cell numbers.
C CALLED FROM:
C COMMUNICATION:
C AUTHOR:        Paul de Jong, NIKHEF-H, 20-02-1990.
C ************************************************************************
C --DATE--:--NAME--:--MODIFICATIONS---------------------------------------
C
C ************************************************************************
C
C PdJ 30-01-91:  The CUCELL table now contains the center of gravity for
C                all cells. In FCAL and RCAL this is equivalent to the
C                center, in the BCAL this might differ a bit (up to 1 cm).
C                For the time being CccXyz gives the center of gravity,
C                but this might change in the future.
C 
C Wuyi Liu 11.02.1998:  shift the BCAL by 6.5 mm towards the RCAL
C                       in the data
C A.Quadt 14.06.1998:   shift the RCAL by 1.1 cm away from 
C                       interaction point in the data
C R.Yoshida 8.02.1998:  Put in all known shifts with respect to the
C                       GAF's from 1996-7.
C
#include "cucell.inc"
#include "partap.inc"
C
#include "cquick.inc"
C
#include "fmckin.inc"
C
#include "zdskey.inc"
C
      INTEGER CellNr
      REAL Cx, Cy, Cz
      INTEGER I, ICNO, IC1, IC2, IPVO, CUR1, CUR2
      REAL RLOW,RHIGH,GETAL1,GETAL2,GETAL3,XR,YR,XT,YT
      LOGICAL FIRST, OK, Error, CExist
      DATA FIRST /.TRUE./
      integer fbr, Modu, Tow, Type, Layer, run_number   
C
      IF (FIRST) THEN
         FIRST=.FALSE.
         ICNO=GETIND(CUCELL,'Nr')

         WRITE(*,*)
         WRITE(*,*) '+-------------------------------+'
         WRITE(*,*) '|                               |'
         WRITE(*,*) '| This is the modified CCCXYZ   |'
         WRITE(*,*) '| F/B/RCAL shifts with respect  |'
         WRITE(*,*) '| to the GAFs are introduced.   |'
         WRITE(*,*) '| 1996-1997 data analysis only  |'
         WRITE(*,*) '|                               |'
         WRITE(*,*) '|     8.2.99  R. Yoshida        |'
         WRITE(*,*) '|                               |'
         WRITE(*,*) '+-------------------------------+'
         WRITE(*,*)
      ENDIF
C
C ..  Does this cell exist?
C
      If (.NOT.CExist(CellNr)) Then
         Error=.TRUE.
         Cx=0.
         Cy=0.
         Cz=0.
         Return
      Endif
C
C
C ..  Fetch its coordinates
C
      Cx = CPosit(CellNr,1)
      Cy = CPosit(CellNr,2)
      Cz = CPosit(CellNr,3)
      Error=.FALSE.
C
C   Now fix the know shifts of the calorimeter with respect to the Gaf's
C      *** This is only valid for 1996-97 data analysis***
C
      If(CouTab(FMCKin).gt.0) return  ! no need to touch MC.
C
C     Get the run number
      run_number   = ZDSKEY_NR1
c
c  Call the routine CellDecode from the EM finder package.
C
c output: FBR   = 0 for FCal, 1 for BCal, 2 for RCal
c         Modu  = Module Number - 1
c         Tow   = Tower Number  - 1
c         Type  = 1-4 EMC, 5 HAC0, 6 HAC1, 7 HAC2
c         Layer = 0 for EMC,HAC0, 1 for HAC1, 2 for HAC2
C
      call CellDecode(CellNr, fbr, Modu, Tow, Type, Layer)
C
C      units are in centimeters. For clarity I have added 1 to the
C      module and tower numbers to get back the usual numbering scheme.
C 
C
      if (fbr.eq.0) then
            Cz = Cz - 0.9
            Cy = Cy + 0.3                   ! FCAL is shifted up by 3 mm up
                                            ! with respect to the GAF's
            if ((Modu+1).lt.12) then
              Cx = Cx - 0.35                ! FCAL south shifted by 3.5 mm
            else if ((Modu+1).eq.12) then
              if ((Tow+1).gt.12) then       ! Mod 12 top is attached to
                 Cx = Cx - 0.35             ! the south.
              else if ((Tow+1).lt.12) then  ! Mod 12 bottom is shifted in
                 Cy = Cy - 0.2              ! y with respect to the rest.
              endif
            endif
      endif
C   BCAL
      if( fbr.eq.1) Cz = Cz - 0.5           ! The Z shift is the compromise
                                            ! value of 5 mm.
c
c   Now RCAL
C
      if (fbr.eq.2) then
            Cz = Cz - 0.9
            Cy = Cy + 0.2                   ! RCAL is shifted up 2 mm
                                            ! with respect to the GAF's
            If (run_number.lt.25596) Cx = Cx - 0.2  ! RCAL shifted 2 mm for 
                                            ! upto run xxxxx
            if ((Modu+1).eq.12) then
              if ((Tow+1).gt.12) then       ! Mod 12 top is shifted in y by 
                 Cy = Cy + 0.1              !  +1 mm
              else if ((Tow+1).lt.12) then  ! Mod 12 bottom is shifted in y
                 Cy = Cy - 0.1              ! by -1 mm
              endif
            endif
      endif           

      RETURN
      END









