c     ====================
      subroutine terminate
c     ====================
c
c --- End of processing termination
c
      implicit none
c
      integer icycle, i
c
c --- Close output file
c
      call hcdir('//histos',' ')
      call hrout(0,icycle,' ')
      call hrend('histos')
      close(50)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)

c --- Delete all histograms and free up memory
c --- Change required for 2022 Linux/Fortran and 2006 Cernlib
      call hdelet(0)
     
      return
      end
