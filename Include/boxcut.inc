C
C 
      box = .false.
      SRTDGOOD = .false.
C
C
      pos_1=-999.0 
      pos_2=-999.0 
      pos_3=-999.0 
C
      cal(1)=e_pos(1)
      cal(2)=e_pos(2)
      cal(3)=e_pos(3)
      cal(4)=elec_ene(1)
C
      srtd(1)=srtd_pos(1)
      srtd(2)=srtd_pos(2)
      srtd(3)=srtd_pos(3)
      srtd(4)=srtd_mip(1)
C
      hes(1)=xhes(1)
      hes(2)=yhes(1)
      hes(3)=-999.0
      hes(4)=ehes(1)
C
      tend(1)=-999.0 
      tend(2)=-999.0 
      tend(3)=-999.0 
      tend(4)=-999.0 
      tend(5)=-999.0 
C
C      CALL best_xyz(pos_1,pos_2,pos_3,cal,srtd,hes,Tend,
C     &     vtx,montecarlo)
C
      whatToUse=1
      CALL best_position(besttheta,pos_1,pos_2,pos_3,cal,srtd,hes
     &     ,Tend,vtx,montecarlo,run,whatToUse,Isys)

C
         if(abs(pos_1).gt.boxcut(1).or.abs(pos_2).gt
     &        .boxcut(2)) then 
           box=.true.
         endif
         if (pos_1.gt.hcut(1).and.pos_1.lt.hcut(2)
     &      .and.abs(pos_2).lt.hcut(5)) then
           box=.false.
         endif
         if (pos_1.gt.hcut(3).and.pos_1.lt.hcut(4)
     &      .and.abs(pos_2).lt.hcut(5)) then
           box=.false.
         endif
C
C
         if (abs(pos_1).gt.outbox(1).or.
     &       abs(pos_2).gt.outbox(2)) then
          box=.false.
         endif
C
