*
*
      INTEGER MaxTrk
      PARAMETER (MaxTrk=10)
      INTEGER ELTfnd, ELTnear, ELTtgid(Maxtrk)
      REAL    ELTmom(4,MaxTrk), ELTsep(MaxTrk), ELTdist(3,MaxTrk)
      REAL    ELTentr(3,MaxTrk), ELTmode(MaxTrk), ELTerr(MaxTrk)
      REAL    ELTang(MaxTrk)
      COMMON /GABFFW/ELTfnd,ELTentr,ELTdist,ELTtgid,ELTmom,ELTsep,
     +               ELTnear,ELTang,ELTmode,ELTerr
*
*
