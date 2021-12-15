c ----------------------------------
      subroutine readdat(choice)
c ----------------------------------

      implicit none

#include "readdata.inc"
      
      integer I
      integer choice
 
      NR_INFILES = 0
      DO I=1,50
	   WRITE(*,*) 'ENTER INPUT NTUPLE NAME : ',I
         READ(CHOICE,90,ERR=95)  FILEN(I)
         IF (LENOCC(FILEN(I)).LT.5) GOTO 95
         NR_INFILES = NR_INFILES + 1
 90      FORMAT(A64)
      ENDDO

 95   continue
	
      WRITE(*,*)
      WRITE(*,*) NR_INFILES,' files read in ...'
      WRITE(*,*)
      WRITE(999,*)
      WRITE(999,*) NR_INFILES,' files read in ...'
      WRITE(999,*)

      return
      end
