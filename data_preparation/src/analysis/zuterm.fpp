C     =================
      SUBROUTINE ZUTERM
C     =================
C
      IMPLICIT NONE
      integer ierr

      WRITE(*,*) ' '
      WRITE(*,*) ' CLOSING NTUPLE '
      WRITE(*,*) ' '

c -----------------------------------
c SRTD termination
c -----------------------------------
      call sreaze(1,0,0,ierr)
      if (ierr.ne.0) write(*,*) 'SREAZE termination error'
c -----------------------------------
c VCTRAK termination
c -----------------------------------
c      call vceaze(-1,1,1,1,1,ierr)
c      if (ierr.ne.0) write(*,*) 'VCEAZE termination error'

      CALL ENDNT

      RETURN
      END

