C     =================
      SUBROUTINE ZUINIT
C     =================
C     ---------------------------------------------------------

      IMPLICIT NONE

#include "sidat95.inc"
c#include "ctrlcm.inc"

      integer ierr,i
      real param(10)

      WRITE(*,*) '                                     '
      WRITE(*,*) '                                     '
      WRITE(*,*) '   +--------------------------------------+'
      WRITE(*,*) '   |                                      |'
      WRITE(*,*) '   |        Ntuple production for         |'
      WRITE(*,*) '   |        F2ISR and FL analysis         |'
      WRITE(*,*) '   |                                      |'
      WRITE(*,*) '   |           Jonathan Scott             |'
      WRITE(*,*) '   |                1999                  |'
      WRITE(*,*) '   +--------------------------------------+'
      WRITE(*,*) '   |    mail to: jpscott@mail.desy.de     |'
      WRITE(*,*) '   |                                      |'
      WRITE(*,*) '   +--------------------------------------+'
      WRITE(*,*) '                                     '
      WRITE(*,*) '                                     '
C

      call booknt

c      call ZESMAF               !! Initialise Magnetic Field

C --- GFLT init ---
c      CALL O1INIT
c      CALL O1ZGID

c      Call GInit
C      Call ZAGeom('ALL', IErr)
C      If (IErr .ne. 0)  Write (*,*) ' ZUINIT: ZAGeom error = ',IErr

c -----------------------------------
c VCTTRAK reconstruction
c -----------------------------------
c      call vceaze(1,0,0,0,1,ierr)
c      if (ierr.ne.0) then
c         write(*,*) 'Error from VCEAZE'
c      end if

      call z_ini_Isles

c      call sreaze(-1,0,0,ierr)
c      if (ierr.ne.0) then
c         write(*,*) 'Error from SREAZE'
c      end if

C --- SIRA95 init ---
      doHESCLU = .TRUE.
      doVCTDCA = .TRUE.
      DCA_eCUT = 25.0

C --- Timing Vtx init ----
      Do I=1,10
         Param(i)=-1.
      EndDo
      Call VtxTimInit(Param)

      RETURN
      END







