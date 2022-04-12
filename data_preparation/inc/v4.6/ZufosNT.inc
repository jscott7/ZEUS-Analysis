C
C v7.0 23/03/00 N.Tuning Add ET
C
      Real               Px,         Py,         Pz,         E,
     +                   Et,
     +                   PxCal,      PyCal,      PzCal,      ECal,
     +                   EtCal,
     +                   PxTrk,      PyTrk,      PzTrk,      ETrk,
     +                   EtTrk,
     +                   PxBSp,      PyBSp,      PzBSp,      EBSp,
     +                   EtBSp
      Real               EminPzNT,   PtNT,       GammaNT,    EtNT

      Common /zRecNT01/  Px,         Py,         Pz,         E,
     +                   Et,
     +                   PxCal(3,2), PyCal(3,2), PzCal(3,2), ECal(3,2),
     +                   EtCal(3,2),
     +                   PxTrk,      PyTrk,      PzTrk,      ETrk,
     +                   EtTrk,
     +                   PxBSp(3),   PyBSp(3),   PzBSp(3),   EBSp(3),
     +                   EtBSp(3),
     +                   EminPzNT,   PtNT,       GammaNT,    EtNT


      Integer   zufos_MaxNT
      Parameter(zufos_MaxNT=250)
      Integer           BSpFlag
      Common /zRecNT02/ BSpFlag(3,zufos_MaxNT)

