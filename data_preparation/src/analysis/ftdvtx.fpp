C
      SUBROUTINE FTDVTX(VTX,MLT)
C     **************************
C ------------------------------------------------------------------------------
C     FTDVTX   Determines the vertex from FTD/RTD SEGs only
C     ======
C
C     Input : Segments from table TFMSEG 
C
C     Output: VTX(1)= 0.
C             VTX(2)= 0.
C             VTX(3)= z_vertex from FTD
C             MLT   = #segments in fit
C
C
C   VERT_FTD=Z_vert= (C^tV^-1C)^-1*(C^tV^-1)*(X_meas)  
C                    \_____________________/
C                                M
C   with         X_fit= (X,Y)_fit= C * VERT_FTD
C
C ------------------------------------------------------------------------------
C   Author: Martin ECKERT PI Bonn
C ------------------------------------------------------------------------------
      IMPLICIT NONE

#include "partap.inc"
#include "tfmseg.inc"


      REAL    VTX(3)
      INTEGER MLT

      INTEGER MAX, NPAR, MAXROW, NCH

      PARAMETER (MAX=50, NPAR=1+2*MAX, MAXROW=1+4*MAX,NCH=50)

      INTEGER I, J, LVMAX, NSEG, COU,
     .        ROWS, COLS, IFAIL, IWORST, IWORST2, ITER,
     .        SID(300), IWIDTH, DENS(NCH), MAXDENS

      REAL    X(MAXROW), W(MAXROW), C(MAXROW,NPAR),
     .        V1(MAXROW),EST(NPAR),
     .        RES(MAXROW), HELP(MAXROW), HLP2(MAXROW), 
     .        CTV(NPAR,MAXROW), CTVC(NPAR,NPAR), CTVCI(NPAR,NPAR),
     .        M(NPAR,MAXROW), RWO(NPAR),R(MAXROW*MAXROW),
     .        XS, YS, ZS, CHI2,
     .        CON(NCH), XWIDTH, ZLIM, ZCLUS, ZDCA(300),
     .        ZFTD1, ZFTD2, ZFTD3, ZRTD
      
      PARAMETER (ZLIM=100.)

      EXTERNAL LVMAX

      LOGICAL FIRST /.TRUE./
C
C --- Code:
C
      IF (FIRST) THEN
        FIRST= .FALSE.
        ZFTD1=  129.
        ZFTD2=  166.
        ZFTD3=  203.
        ZRTD = -129. 

        CALL HBOOK1(111111,'zdh',NCH,-ZLIM,ZLIM,0.)
      ENDIF  

      ITER= 0
      NSEG= 0
      COU = 0
      MLT = 0

      CALL VZERO(C,MAXROW*NPAR)
      CALL VZERO(VTX,3)

      CALL HRESET(111111,' ')

C --- Loop over FTD1/RTD SEGs and select good guys:
      
      DO 50 I=1, COUTAB(TFMSEG)
        CALL FETTAB(TFMSEG,ID,I)

C -     SEGments in FTD1,2,RTD:
        IF ( ABS(TFMSEG_ZS-ZFTD1) .LT. 5. .OR.
     .       ABS(TFMSEG_ZS-ZFTD2) .LT. 5. .OR.
     .       ABS(TFMSEG_ZS-ZRTD ) .LT. 5.     ) THEN

C ---   Check if SEG can be used for vertexing (z position of closest approach)
          ZS= -(TFMSEG_xs*TFMSEG_xsp + TFMSEG_ys*TFMSEG_ysp)
     .        /(TFMSEG_xsp**2        + TFMSEG_ysp**2)
     .        + TFMSEG_zs 

C ---     Its distance to beam pipe:
          XS= TFMSEG_xs+TFMSEG_xsp*(ZS-TFMSEG_zs)
          YS= TFMSEG_ys+TFMSEG_ysp*(ZS-TFMSEG_zs)         

          IF (SQRT(XS**2+YS**2) .LT.  10. .AND.
     .        ABS(ZS)           .LT. ZLIM      ) THEN

C -         It's OK, fill in histo & remember:
            CALL HF1(111111,ZS,1.)

            COU      = COU+1
            SID (COU)= TFMSEG_ID 
            ZDCA(COU)= ZS
          ENDIF

        ENDIF  
 50   ENDDO

C
C --- Examine histo for clusters:
C

      CALL VZERO(DENS,NCH)
      CALL HUNPAK(111111,CON,' ',1)

      XWIDTH= 20.
      IWIDTH= (XWIDTH*NCH)/(2*ZLIM)
      DO I= IWIDTH/2+1, NCH-IWIDTH/2
        DO J= -(IWIDTH)/2, IWIDTH/2
          DENS(I)= DENS(I)+CON(I+J) 
        ENDDO          
      ENDDO  

      MAXDENS= LVMAX(DENS,NCH-IWIDTH)
      MAXDENS= MAXDENS
      CALL HIX(111111,MAXDENS,ZCLUS)

c      WRITE(*,*)'Cluster: ',MAXDENS,ZCLUS
c      CALL HPRINT(111111)

C --- Now handle selected SEGs:

      DO 100 J= 1, COU

C -     Distance of closest approach :
c        IF (ABS(ZDCA(J)-ZCLUS).LT.25.) THEN
        IF (ABS(ZDCA(J)-ZCLUS).LT.15.) THEN

          TFMSEG_ID= SID(J)
          CALL GETTAB(TFMSEG)
          
          NSEG= NSEG+1
          ROWS= 1+4*NSEG
          COLS= 1+2*NSEG


C --- Fill Vector X :
c     ((X,Y) measurements of SEGs @z_seg) 

C ---     nominal vertex:
          X(1)     = 0.

          X(ROWS-3)= TFMSEG_xs
          X(ROWS-2)= TFMSEG_ys
          X(ROWS-1)= TFMSEG_xsp
          X(ROWS  )= TFMSEG_ysp
          

C --- Fill Matrix C :
c     (describes how to get the (X,Y)_extr from the 3 vertex infos &
c      refitted slopes from original slopes)

          C(1,1)=  1.

          C(ROWS-3,1)=  -TFMSEG_xsp
          C(ROWS-2,1)=  -TFMSEG_ysp
          
          C(ROWS-3,COLS-1)= TFMSEG_zs
          C(ROWS-3,COLS  )= 0.

          C(ROWS-2,COLS-1)= 0.
          C(ROWS-2,COLS  )= TFMSEG_zs

          C(ROWS-1,COLS-1)= 1.
          C(ROWS-1,COLS  )= 0.

          C(ROWS  ,COLS-1)= 0.
          C(ROWS  ,COLS  )= 1.

C --- Inverse Covariance Matrix (diagonal elements only) V^-1

          V1(1)     = 1.
          V1(ROWS-3)= 1./TFMSEG_cov(1 )
          V1(ROWS-2)= 1./TFMSEG_cov(8 )
          V1(ROWS-1)= 1./TFMSEG_cov(5 )
          V1(ROWS  )= 1./TFMSEG_cov(10)

        ENDIF

        IF (NSEG .EQ. MAX) THEN
          WRITE(*,*)' @FTDVERTEX: Too many SEGs!', COU
          GOTO 101
        ENDIF  


 100  ENDDO  
 101  CONTINUE

      IF (NSEG.LT.3)           GOTO 999

C --- Iteration label:
 123  CONTINUE

      CALL VZERO(CTV  ,NPAR*MAXROW)
      CALL VZERO(CTVC ,NPAR*NPAR)
      CALL VZERO(CTVCI,NPAR*NPAR) 
      CALL VZERO(RWO,NPAR)
      
      CALL VZERO(M,NPAR*MAXROW)  
      CALL VZERO(EST,NPAR)

      CALL VZERO(HELP,MAXROW)
      CALL VZERO(HLP2,MAX)
      CALL VZERO(RES,MAXROW)


C --- Update values (z_vert only):

      CALL UCOPY(W,X,1)
      
C --- 1. C^t * V^-1 = CTV

      DO I= 1, COLS
        DO J= 1, ROWS
          CTV(I,J)= C(J,I)*V1(J) 
        ENDDO
      ENDDO
  
C --- 2. CTV * C   = CTVC

      CALL RMMLT(COLS,ROWS,COLS,
     .           CTV,CTV(1,2),CTV(2,1),
     .           C,C(1,2),C(2,1),
     .           CTVC,CTVC(1,2),CTVC(2,1),
     .           R)

C --- 3. CTVC^(-1) = CTVCI 


      DO I= 1, COLS
        CALL UCOPY(CTVC(1,I),CTVCI(1,I),COLS)        
      ENDDO  

      CALL RINV(COLS,CTVCI,NPAR,RWO,IFAIL)
 
      IF (IFAIL .NE. 0) THEN
        WRITE(*,'(A,I3)') 'Matrix inversion failed!, NSEG=',NSEG
        GOTO 999
      ENDIF

C --- 4. CTVCI * CTV -> M

      CALL RMMLT(COLS,COLS,ROWS,
     .           CTVCI,CTVCI(1,2),CTVCI(2,1),
     .           CTV,CTV(1,2),CTV(2,1),
     .           M,M(1,2),M(2,1),
     .           R)


C --- 5. The fit itself :
c     (by multiplication of Estimate matrix M with measurement vector X)
c      EST= M * X


      CALL RMMLT(COLS,ROWS,1,
     .           M,M(1,2),M(2,1),
     .           X,X(1),X(2),
     .           EST,EST(1),EST(2),
     .           R)

C --- 6. Chisquare of fit:
c        C*EST = (X,Y)_extr= W()

      CALL RMMLT(ROWS,COLS,1,
     .           C,C(1,2),C(2,1),
     .           EST,EST(1),EST(2),
     .           W,W(1),W(2),
     .           R)

c -   RES= (X,Y)_meas - (X,Y)_fit = X()-W()

      CALL RVSUB(ROWS,X(1),X(2),W(1),W(2),RES(1),RES(2))

c -   CHI2= RES^t * V^-1 * RES   

c      CALL RMMLT(ROWS,ROWS,1,
c     .           V,V(1,2),V(2,1),
c     .           RES,RES(1),RES(2),
c     .           HELP,HELP(1),HELP(2))


c      CALL RMMLT(1,ROWS,1,
c     .           RES,RES(2),RES(1),
c     .           HELP,HELP(1),HELP(2),
c     .           vert_res,vert_res,vert_res)


C -   CHI2 old fashioned      

c     VTX residuals:
c      CHI2= RES(1)**2*V1(1)+RES(2)**2*V1(2)+RES(3)**2*V1(3)
      CHI2= RES(1)**2*V1(1)

c     SEG residuals:
      DO I= 1, (ROWS-1)/4
        J=4*I+1
        HELP(J-4)= RES(J-3)**2*V1(J-3)
        HELP(J-3)= RES(J-2)**2*V1(J-2)
        HELP(J-2)= RES(J-1)**2*V1(J-1)
        HELP(J-1)= RES(J  )**2*V1(J  )

        HLP2(I)  = HELP(J-4)+HELP(J-3)+HELP(J-2)+HELP(J-1)
        CHI2= CHI2 + HLP2(I)
      ENDDO  

C-    Chi2 / NDF:
c      vert_res= CHI2/(ROWS-COLS)

C --- Output:

c      vert_ftde= SQRT(CTVCI(1,1))


C --- Check fit result, good enough?

C -   Too bad -> find worst SEG and remove

        ITER= ITER+1

        IWORST2= LVMAX(HLP2,NSEG)

c        IF (vert_res       .GT. 20. .OR.
c     .      HLP2(IWORST2)  .GT. 25.      ) THEN 

C-    Chi2 / NDF:

        IF (CHI2/(ROWS-COLS) .GT. 20. .OR.
     .      HLP2(IWORST2)    .GT. 25.      ) THEN 

C -     Shift array entries:

        IWORST= IWORST2*4-2

        IF (IWORST2.LT.MAX) THEN
          
          DO I= IWORST+4, ROWS
            X (I-4)  = X (I)
            X (I)    = 0. 
            V1(I-4)  = V1(I)
            V1(I)    = 0.
            C (I-4,1)= C (I,1)
            C (I,  1)= 0.

            DO J= 4, COLS
              C(I-4,J-2)= C(I,J)
              C(I  ,J  )= 0.
            ENDDO  
          ENDDO  

        ENDIF

c -     Decrement counters:
        ROWS= ROWS-4
        COLS= COLS-2
        NSEG= NSEG-1

C -     Next Iteration
        IF (NSEG.GT.2) GOTO 123

      ENDIF  


C --- Output values:

      VTX(1)= 0.
      VTX(2)= 0.
      VTX(3)= EST(1)
      MLT   = NSEG

C ---
 999  RETURN
      END

******************************************************************************
C
