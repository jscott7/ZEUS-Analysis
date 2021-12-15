C
C
C------------------------------------------------------------------
C
      SUBROUTINE  PDFINIT
c     --------------------
      IMPLICIT NONE 

      REAL q2,x,f2data,f2mc,wei,rewrap
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VAL(20)
      integer g
      logical first
      save first
      data first /.true./

      if (first) then
         WRITE (*,*) ' ****************************'
         WRITE (*,*) ' '
         WRITE (*,*) '         PDF INIT '
         WRITE (*,*) ' '
         WRITE (*,*) ' ****************************'
         WRITE (999,*) ' ****************************'
         WRITE (999,*) ' '
         WRITE (999,*) '         PDF INIT '
         WRITE (999,*) ' '
         WRITE (999,*) ' ****************************'


c Fill common blocks
         PARM(1) = 'Init0'
         VAL(1)  = 0.D0
         CALL PDFSET(PARM,VAL)

c Particle type 1 = nucleons
         PARM(1) = 'Nptype'
         VAL(1)  = 1

c Author group 3 = MRS
         PARM(2) = 'Ngroup'
         VAL(2)  = 3

c Structure Function set 44 = MRSA Low Q^2
         PARM(3) = 'Nset'
         VAL(3)  = 44

         CALL PDFSET(PARM,VAL)
         first=.false.
      endif

      return
      end

C
C
C------------------------------------------------------------------
C
C     ------------------------
      SUBROUTINE CF2 (x,q2,F2)
c     ------------------------
      IMPLICIT NONE
      REAL*4  XPQ(-6:6),F2,CHARGE(6),x,q2
      DATA    CHARGE /-0.33333,0.66667,-0.33333,
     &     0.66667,-0.33333,0.66667/
      INTEGER I
      double precision dx,dq,upv,dnv,ups,dns,sts,chs,bos,tos,glu
      dx=dble(x)
      dq=dble(sqrt(max(0.,q2)))

c THIS WAS A BUG!
c      dq2=q2
      call structm(dx,dq,upv,dnv,ups,dns,sts,chs,bos,tos,glu)
      xpq( 1)=dnv+dns
      xpq(-1)=dns
      xpq( 2)=upv+ups
      xpq(-2)=ups
      xpq( 3)=sts
      xpq(-3)=sts
      xpq( 4)=chs
      xpq(-4)=chs
      xpq( 5)=bos
      xpq(-5)=bos
      xpq( 6)=tos
      xpq(-6)=tos
      F2 = 0.0
      DO I = 1,6
         F2 = F2 + CHARGE(I)**2*(XPQ(I)+XPQ(-I))
      ENDDO
      RETURN
      END
C


