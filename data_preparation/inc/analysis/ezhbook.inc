      INTEGER NWDS_HBOOK
      PARAMETER (NWDS_HBOOK=1200000)
      REAL          HMEM
      COMMON /PAWC/ HMEM(NWDS_HBOOK)

      CHARACTER*132 NtName
      INTEGER       ntLun, NT_ID
      DATA          ntLun / 10 /
      DATA          NT_ID / 10 /


