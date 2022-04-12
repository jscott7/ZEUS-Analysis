

*-->  CAL geomerty variables
      Integer     get_fbr, get_mod, get_tow, get_typ
      Integer     cellincode

C     These are intrinsic functions:
      integer     ishft , iand

*-->  define some CAL Information functions (F=0,B=1,R=2)
      get_fbr(Cell) = ISHFT(Cell,-14)
      get_Mod(Cell) = ISHFT(IAND(Cell,x'3E00'),-9)
      get_Tow(Cell) = ISHFT(IAND(Cell,x'1F0' ),-4)
      get_typ(Cell) = ISHFT(IAND(Cell,x'E'   ),-1)

      cellincode(CAL,Module,Tower,cType) = 
     &     16384*CAL+(Module-1)*512+(Tower-1)*16+cType*2
*
