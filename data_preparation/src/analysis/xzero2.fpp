C     =========================================
      Real Function Xzero2(zvtx,theta,phi,year)
C     =========================================
c
      Implicit None
c
      Real zvtx,theta,phi,dz,dtheta,dphi,zMax
      integer nPhi, nTheta, nZ, nYear, year
      parameter (nPhi=72, nTheta=720, nZ=3, zMax = 45.0, nYear=3)
      Dimension dead94(nZ,nTheta,nPhi),dead95(nZ,nTheta,nPhi)
      Dimension dead96(nZ,nTheta,nPhi)
      Real dead94,dead95,dead96,distz,disttheta,distphi,dist
      Real pi,Zbin(2),thetabin(2),phibin(2),denom
      Integer i,j,k,Icall,IZ(2),Itheta(2),Iphi(2)
      character FileName(nZ*nYear)*14
      Data Icall/0/,pi/3.14159/
      integer newlun
      logical ifex
      DATA FILENAME/'deadm45_94.dat','dead0_94.dat','deadp45_94.dat',
     +     'deadm45_95.dat','dead0_95.dat','deadp45_95.dat',
     +     'deadm45_96.dat','dead0_96.dat','deadp45_96.dat'/
c
      Icall = Icall + 1
c
c  Initialize
      If (Icall.eq.1) then
         PRINT *,' '
         PRINT *,'This is a new Xzero. It can use either 94, 95'//
     +        ' or 96 X0-map!'
         PRINT *,' '
         DO newlun=21,199
            INQUIRE (UNIT=newlun, IOSTAT=I, ERR=7001,OPENED=IFEX)
            IF(.NOT.IFEX) goto 7101
         ENDDO
c     opening error
 7001    PRINT *,'error',I,' in Xzero2; couldn''t find free LUN'
         STOP
c     whow! got a lun to proceed ...
 7101    CONTINUE
         DO K=1,nZ
            inquire( file=FILENAME(K),exist=IFEX )
            IF(IFEX) THEN
               open(newlun,file=FILENAME(K),status='old')
               Read(newlun,800) ((dead94(K,i,j),i=1,nTheta),j=1,nPhi)
               close(newlun)
            ELSE
               PRINT *,'couldn''t find file '//FILENAME(K)
               STOP
            ENDIF
            inquire( file=FILENAME(K+nYear),exist=IFEX )
            IF(IFEX) THEN
               open(newlun,file=FILENAME(K),status='old')
               Read(newlun,800) ((dead95(K,i,j),i=1,nTheta),j=1,nPhi)
               close(newlun)
            ELSE
               PRINT *,'couldn''t find file '//FILENAME(K+nYear)
               STOP
            ENDIF
            inquire( file=FILENAME(K+2*nYear),exist=IFEX )
            IF(IFEX) THEN
               open(newlun,file=FILENAME(K),status='old')
               Read(newlun,800) ((dead96(K,i,j),i=1,nTheta),j=1,nPhi)
               close(newlun)
            ELSE
               PRINT *,'couldn''t find file '//FILENAME(K+2*nYear)
               STOP
            ENDIF
         ENDDO
c
         dtheta = pi/float(nTheta)
         dphi   = 2.*pi/float(nPhi)
         dz = 2.*zMax/float(nZ-1)
      Endif
c
 800  Format(10(1x,F6.3))
c
c  Get indices
c
      IZ(1) = 1 + INT((zvtx+zMax)/dz)
      IZ(1) = max(IZ(1),1)
      IZ(1) = min(IZ(1),2)
      IZ(2) = IZ(1) + 1
      Zbin(1)  = -zMax + float(IZ(1)-1)*dz
      Zbin(2)  = Zbin(1) + dz
c
      Itheta(1) = INT( theta/dtheta ) + 1
      Itheta(1) = max(Itheta(1),1)
      Itheta(1) = min(Itheta(1),nTheta-1)
      Itheta(2) = Itheta(1) + 1
      thetabin(1)  = (Itheta(1)-1)*dtheta
      thetabin(2)  = thetabin(1) + dtheta
c
      Iphi(1) = Mod(INT(phi/dphi+float(2*nPhi)),nPhi) + 1
      Iphi(1) = max(Iphi(1),1)
      Iphi(1) = min(Iphi(1),nPhi-1)
      Iphi(2) = Iphi(1) + 1
      phibin(1)  = (Iphi(1)-1)*dphi
      phibin(2)  = phibin(1) + dphi
c
      Xzero2 = 0.
      denom = 0.
c
      Do i=1,2
         distz = (zvtx - Zbin(i))/dz
         Do j=1,2
            disttheta = (theta - thetabin(j))/dtheta
            Do k=1,2
               distphi = (phi - phibin(k))/dphi
               dist = sqrt(distz**2 + disttheta**2 + distphi**2)
               IF(year.EQ.1994) THEN
                  Xzero2 = Xzero2 + dead94(IZ(i),Itheta(j),Iphi(k))/
     +                 dist**2
               ELSEIF(year.EQ.1995) THEN
                  Xzero2 = Xzero2 + dead95(IZ(i),Itheta(j),Iphi(k))/
     +                 dist**2
               ELSEIF(year.EQ.1996.OR.year.EQ.1997) THEN
                  Xzero2 = Xzero2 + dead96(IZ(i),Itheta(j),Iphi(k))/
     +                 dist**2
               ELSE
                  PRINT *,'WARNING from Xzero: wrong year selected: ',
     +                 year
                  Xzero2 = 0.
                  RETURN
               ENDIF
               denom = denom + 1/dist**2
            Enddo
         Enddo
      Enddo
c
      Xzero2 = Xzero2/denom
      End
