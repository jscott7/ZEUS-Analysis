c -------------------------------------
      subroutine f2wt97
c -------------------------------------

	implicit none

#include"local.inc"
#include"common97.inc"

	logical first
      data first/.true./
      save first
      logical first2
      data first2/.true./
      save first2

	fallm = 0.
	f2mrsa = 0.
	f2weight = 0.

c ----------------------------------
c Calculate f2 values 
c ----------------------------------
      if (first) then
         write(*,*) '*************************'
         write(*,*) '*** calculate f2 allm ***'
         write(*,*) '*************************'
         write(999,*) '*************************'
         write(999,*) '*** calculate f2 allm ***'
         write(999,*) '*************************'

         first = .false.
      end if
      call f2allm(x_tru, q2_tru, fallm)

      if (first2) then
         write(*,*) '*************************'
         write(*,*) '*** calculate f2 MRSA ***'
         write(*,*) '*************************'
         write(999,*) '*************************'
         write(999,*) '*** calculate f2 MRSA ***'
         write(999,*) '*************************'

         first2 = .false.
      end if
      call CF2(x_tru, q2_tru, f2mrsa)

c -----------------------------------------
c Weight MC 
c -----------------------------------------
      f2weight = fallm / f2mrsa

	return
	end
