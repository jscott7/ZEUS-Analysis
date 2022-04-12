 
      INTEGER           MAXbm, IGMUON
      PARAMETER        (MAXbm = 63)
      COMMON / MUOTUP / Nbmtk,
     &                  nRefAX(MAXbm),  Xseg (MAXbm),
     &                  Yseg  (MAXbm),  Zseg (MAXbm),
     &                  XDCos (MAXbm),  YDCos(MAXbm),
     &                  ZDCos (MAXbm),  InvP (MAXbm),
     &                  BMChi2(MAXbm),  IGMUON
 
      INTEGER           Nbmtk, nRefAX
      REAL              Xseg,  Yseg,  Zseg,
     &                  XDCos, YDCos, ZDCos,
     &                  InvP,  BMChi2
 
      CHARACTER*(500)   MUOTUP_form(2)
      DATA              MUOTUP_form /
     &                  'MUOTUP',
     &'Nbmtk:I::[0,63],nRefAX(Nbmtk):I::[1,9],Xseg(Nbmtk):R,
     &Yseg(Nbmtk):R,Zseg(Nbmtk):R,XDCos(Nbmtk):R,YDCos(Nbmtk):R,
     &ZDCos(Nbmtk):R,InvP(Nbmtk):R,BMChi2(Nbmtk):R,IGMUON:I' /
 
 
 
 
 
