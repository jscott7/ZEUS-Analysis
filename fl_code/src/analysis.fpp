c ----------------------------
      program main
c ----------------------------
      implicit none

#include "readdata.inc"
#include "local.inc"
	
      INTEGER I,CHOICE,num
      integer histofill
      real x1(50)
      real f21(50)

C ------------------------------
C --- INPUT NTUPLE FILENAMES ---
C ------------------------------
      WRITE(*,*) 'READ FROM :  '
      WRITE(*,*) '------------------------------------------'
      WRITE(*,*) '            1 - Do 96 FL/ISR'
      WRITE(*,*) '            2 - Do 97 FL/ISR'
      WRITE(*,*) '            3 - Do all FL/ISR'
      WRITE(*,*) '------------------------------------------'
      READ(*,*) CHOICE

      num = 1
      i = 1
      if (choice.eq.3) num=26
      do i = 1,num
c ------------------------------------------
c Initialisation
c ------------------------------------------
        call isr_init(choice,i)
	
        call crumpet

c ------------------------------------------
c Run over all relevant Ntuples.
c ------------------------------------------
	if (choice.eq.1) then 
	   call do_all96(i)
	else if (choice.eq.2) then
	   call do_all97(i)
	else if (choice.eq.3) then
	   call do_all(i)
	end if

c --------------------------------
c Write out bins
c --------------------------------
	write(*,*) '************************'
	write(*,*) '*** Writing out Bins ***'
	write(*,*) '************************'
	write(999,*) '************************'
	write(999,*) '*** Writing out Bins ***'
	write(999,*) '************************'

	write(81,*) binbgd
	write(82,*) binmc
	write(83,*) bindat
	write(84,*) f2isr
	write(85,*) errors
		
	call terminate
	end do

	END

c Math errors handling routine
c      SUBROUTINE MATHERRQQ( name, length, info, retcode)
c	USE DFLIB
c	 INTEGER(2) length, retcode
c	 CHARACTER(length) name
c	RECORD /MTH$E_INFO/ info
c	  PRINT *, "Entered MATHERRQQ"
c	PRINT *, "Failing function is: ", name
c	 PRINT *, "Error type is: ", info.errcode
c	 IF ((info.ftype == TY$REAL4 ).OR.(info.ftype == TY$REAL8)) THEN
c      PRINT *, "Type: REAL"
c     PRINT *, "Enter the desired function result: "
c    READ(*,*) info.r8res
c       retcode = 1
c	 ELSE IF((info.ftype==TY$CMPLX8 ).OR.
c    &(info.ftype==TY$CMPLX16)) THEN
c       PRINT *, "Type: COMPLEX"
c       PRINT *, "Enter the desired function result: "
c       READ(*,*) info.c16res
c       retcode = 1
c	 END IF
c	END
