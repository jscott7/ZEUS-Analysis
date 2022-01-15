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
     
      return
      end
