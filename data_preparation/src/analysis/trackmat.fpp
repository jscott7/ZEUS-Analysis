c     ========================================================
      subroutine trackmat(calpos,rcut,track_p,r_dist,
     &                    trkpos,trk_iso1,trk_iso2,index,ierr)
c     ========================================================
c
      implicit none
c     -------------------------------------------
c
c >>> Inputs
c
      real calpos(3)
c     ==============
c --- Contains the calorimeter position
c
      real rcut
c     =========
c --- Cut on separation between extrapolated track and
c     cal position (in cm). r defines a sphere centred
c     on the cal object.
c
c     -------------------------------------------
c
c >>> Outputs
c
      real track_p
c     ============
c --- Momentum of closest track to given position
c
C     dist (in cm)
C     trkpos(3) x,y,z pos of track extrapolation
C     trk_iso1, trk_iso2
C
C     integer index -> VCTRHL track ID
C

      integer ierr
c     ============
c --- Status return
c
c     0 -> All OK
c     1 -> No track within rcut cm of position
c
c     -------------------------------------------
C
#include "partap.inc"
#include "vcatcal.inc"
#include "vctrhl.inc"
c
      integer i,j,k,isel
      logical first,debug
      save first
      data first/.true./
      data debug/.false./


c ----------------------------
c --- Loose cut definition ---
c ----------------------------
      real loosecut
      parameter (loosecut = 50.)
c
      real r,r2,l(3),a(3),n(2,2),lambda,xpos(2,2),dist,xpos3(3),
     +     spos3(3),l3(3),p,q,nenner,dist2,r_sum, r_min,
     +     trk_iso1,trk_iso2,r_dist,trk_x,trk_y,trk_z,trk_p,
     +     trkpos(3),bcal

      integer index, n_trk, n_cand
c
      IF (debug) WRITE(*,*)
      IF (debug) WRITE(*,*) 'THIS IS TRACKMAT ...'
      IF (debug) WRITE(*,*) '===================='

C -------------------------------
C --- Init tracking variables ---
C -------------------------------
      trk_x     = -9999.
      trk_y     = -9999.
      trk_z     = -9999.
      trk_p     = -9999.
      bcal      = -9999.
      n_cand    = -9999.
      trk_iso1  = -9999
      trk_iso2  = -9999.
      r_dist    = -9999.
      trkpos(1) = -9999.
      trkpos(2) = -9999.
      trkpos(3) = -9999.


      if (first) then
         first = .false.
         call cresel(VCATCAL,isel,'trkmatch')
      endif
c
      if (coutab(VCATCAL).lt.1) then
         ierr = 1
         return
      endif
c
      call clesel(VCATCAL,isel)

c ----------------------------------------------------------
c --- First do a loose selection in the area of interest ---
c ----------------------------------------------------------
      IF (debug) write(*,*) 'CAL_POS IS :',calpos(1),calpos(2),calpos(3)
      IF (debug) WRITE(*,*) 'do preselection ...'

      do i=1,coutab(VCATCAL)
         call fettab(VCATCAL,ID,i)
         r = sqrt((VCATCAL_X-calpos(1))**2 +
     &            (VCATCAL_Y-calpos(2))**2 +
     &            (VCATCAL_Z-calpos(3))**2)
         if (r.lt.loosecut) call inssel(VCATCAL,isel)
         IF (debug) WRITE(*,*) VCATCAL_X,VCATCAL_Y,VCATCAL_Z,r
      enddo

      n_cand = float(cousel(VCATCAL,isel))

c ------------------------------------------------------------
c --- Need at least one track passing loose selection cuts ---
c ------------------------------------------------------------
      if (cousel(VCATCAL,isel).lt.1) then
         ierr = 1
         return
      endif
c
      if (cousel(VCATCAL,isel).gt.1) then
c
C --------------------------------------------------------------------
c --- Now loop over selected tracks doing a straight line          ---
c --- extrapolation to the point of closest approach to the        ---
c --- cal position. If more than 1 track is selected then take the ---
c --- closest one.                                                 ---
c --------------------------------------------------------------------
         dist  = 1000.
         index = -1
         do i=1,cousel(VCATCAL,isel)
            call fettab(VCATCAL,isel,i)

c ----------------------------
c --- Extrapolation vector ---
c ----------------------------
            l(1) = VCATCAL_PX
            l(2) = VCATCAL_PY
            l(3) = VCATCAL_PZ
c
c --- Point to extrapolate from
c
            a(1) = VCATCAL_X
            a(2) = VCATCAL_Y
            a(3) = VCATCAL_Z
c
c --- extrapolate to CAL ---
c
            IF (calpos(3).lt.-148.) THEN    !!! RCAL
               IF (l(3).NE.0.0) THEN
                  lambda = (calpos(3) - a(3)) / l(3)
               ELSE
                  lambda = -9999.
               ENDIF
            ELSE                            !!! BCAL
               nenner = (l(1)**2 + l(2)**2)
               IF (nenner.gt.0.) THEN
                  p = 2.*(a(1)*l(1) + a(2)*l(2))/nenner
                  q = (a(1)**2+a(2)**2-calpos(1)**2-calpos(2)**2)
     &                / nenner
                  IF (0.25*p**2-q.gt.0.0) THEN
                    lambda = -0.5*p + sqrt(0.25*p**2 - q)
                  ELSE
                    lambda = -9999.
                  ENDIF
               ELSE
                  lambda = -9999.
               ENDIF
            ENDIF

            do k=1,3
               xpos3(k) = a(k)+(lambda*l(k))
            enddo
c
c --- xpos3 contains the extrapolated position of closest approach
c     of the track to the cal position
c
            r = sqrt( (xpos3(1)-calpos(1))**2+
     +           (xpos3(2)-calpos(2))**2+
     +           (xpos3(3)-calpos(3))**2 )
c
            if (r.lt.dist) then
               dist  = r
               index = i
            endif
c
         enddo
c
      else
         index = 1
      endif

C --------------------------------------------------------
c --- Get parameters for best matching track -------------
c --- Now do the exercise again for the closest track  ---
c --- (lots of code but should be fast)                ---
c --------------------------------------------------------

      call fettab(VCATCAL,isel,index)
c
c --- Extrapolation vector
c

      l(1) = VCATCAL_PX
      l(2) = VCATCAL_PY
      l(3) = VCATCAL_PZ
c
c --- Point to extrapolate from
c
      a(1) = VCATCAL_X
      a(2) = VCATCAL_Y
      a(3) = VCATCAL_Z
c
c --- Extrapolate to CAL
c
            IF (calpos(3).lt.-148.) THEN    !!! RCAL
               bcal = 0.0
               IF (l(3).NE.0.0) THEN
                  lambda = (calpos(3) - a(3)) / l(3)
               ELSE
                  lambda = -9999.
               ENDIF
            ELSE                            !!! BCAL
               bcal   = 1.0
               nenner = (l(1)**2 + l(2)**2)
               IF (nenner.gt.0.) THEN
                  p = 2.*(a(1)*l(1) + a(2)*l(2))/nenner
                  q = (a(1)**2+a(2)**2-calpos(1)**2-calpos(2)**2)
     &                / nenner
                  IF (0.25*p**2-q.gt.0.0) THEN
                    lambda = -0.5*p + sqrt(0.25*p**2 - q)
                  ELSE
                    lambda = -9999.
                  ENDIF
               ELSE
                  lambda = -9999.
               ENDIF
            ENDIF

      do k=1,3
         xpos3(k) = a(k)+(lambda*l(k))
      enddo

c ------------------------------------------------------------
c --- Now really do have the track we want and its gubbins ---
c ------------------------------------------------------------
      if (debug) write(*,*) 'tracking pos: ',xpos3(1),xpos3(2),xpos3(3)
      r = sqrt( (xpos3(1)-calpos(1))**2+
     +          (xpos3(2)-calpos(2))**2+
     +          (xpos3(3)-calpos(3))**2 )
c
      r_dist = r

      if (r.lt.rcut) then
         index   = VCATCAL_ID
         call fettab(VCTRHL,ID,index)
         track_p = VCTRHL_pgevc
         ierr    = 0

C ------------------------------
C ---    Determine isolation ---
C ------------------------------
         n_trk = 0
         r_sum = 0.0
         r_min = 9999.
         do i=1,coutab(VCATCAL)
            call fettab(VCATCAL,ID,i)
            r = sqrt((VCATCAL_X-calpos(1))**2 +
     &               (VCATCAL_Y-calpos(2))**2 +
     &               (VCATCAL_Z-calpos(3))**2)
            IF (i.NE.index) THEN
               r_sum = r_sum + r
               n_trk = n_trk + 1
               IF (r.LT.r_min) r_min = r
            ENDIF
         enddo
      else
         ierr = 1
      endif

C -------------------------
C --- output parameters ---
C -------------------------

      if (ierr.eq.0) THEN
        if (debug) write(*,*) 'track momentum : ',track_p,
     &             ' track end point: ',trk_x,trk_y,trk_z,
     &             ' cal position : ',(calpos(i),i=1,3)

        trk_x     = xpos3(1)
        trk_y     = xpos3(2)
        trk_z     = xpos3(3)
        trk_p     = track_p
        trkpos(1) = trk_x
        trkpos(2) = trk_y
        trkpos(3) = trk_z
        IF (n_trk.gt.0) THEN
          trk_iso1 = r_sum / float(n_trk)
        ELSE
          trk_iso1 = 9999.
        ENDIF
        trk_iso2 = r_min
      else
        if (debug) write(*,*) 'no track matched'
      endif
c
      return
      end

