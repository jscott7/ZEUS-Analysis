*.............................................................................
*  (cell/cone)Island Information.............
*  This KEEP can be used with both cellIsland and coneIsland routines.
*  When making coneIslands: cellIslands are made first and then 
*  recombined into coneIslands based on the criterion of clossest 
*  approach in theta-phi space.
*  For more info look at http://www-zeus.desy.de/~briskin/ZEUS/zufos.html
* 
*  cellIslands ==> routine z_cIsles
*  coneIslands ==> routine Isles      
*.............................................................................

*...
*  Nr_cMax    <==> Maximum number of cells in the event
*  nI_Max     <==> Maximum number of (cell)/(cone)Islands in the event
*  Nr_cIslMax <==> Maximum number of cells in the (cell)/(cone)Island
*...
      Integer    Nr_cMax      , nI_Max     , Nr_cIslMax
      Parameter( Nr_cMax=1000 , nI_Max=200 , Nr_cIslMax=250)
*...
*  Copy of Caltru Information........
*  zCells            <==> Number    of cells to be considered for clustering.
*  zPnrl( 1:Nr_cMax) <==> PoserNr   of the cell
*  zEl(   1:Nr_cMax) <==> Energy    of the cell
*  zImbl( 1:Nr_cMax) <==> Imbalance of the cell
*  zT1(   1:Nr_cMax) <==> time(1)   of the cell
*  zT2(   1:Nr_cMax) <==> time(2)   of the cell
*  zID(   1:Nr_cMax) <==> Row number in table Caltru(i.e. Caltru_ID)
*...
      Integer        zCells,zPnrl,zID
      Real           zEl   ,zImbl,zT1,zT2
      Common/zIsl01/ zCells,zPnrl(Nr_cMax),
     +                      zEl(  Nr_cMax),
     +                      zImbl(Nr_cMax),
     +                      zT1(  Nr_cMax),
     +                      zT2(  Nr_cMax),
     +                      zID(  Nr_cMax)

*...
*  nIsl              <==> Number of (cell/cone)Islands after clustering
*  NrcIsl(  1:nIsl)  <==> Number of cells in (cell/cone)Island
*  eIsl(    1:nIsl)  <==> Energy of (cell/cone)Island
*  xIsl(    1:nIsl)  <==> X Center-of-Gravity (COG) using log weighting
*  yIsl(    1:nIsl)  <==> Y Center-of-Gravity (COG) using log weighting
*  zIsl(    1:nIsl)  <==> Z Center-of-Gravity (COG) using log weighting
*  rIsl(    1:nIsl)  <==> Maximum Radius of (cell/cone)Island              
*  emcEIsl( 1:nIsl)  <==> EMC section Energy of (cell/cone)Island
*  xI1(     1:nIsl)  <==> X Center-of-Gravity (COG) using liniar weighting
*  yI1(     1:nIsl)  <==> Y Center-of-Gravity (COG) using liniar weighting
*  zI1(     1:nIsl)  <==> Z Center-of-Gravity (COG) using liniar weighting
*
*  Internal information for Isles routine.......
*  IslTyp           <==> This variable makes sence only 
*                        when using Isles routine with
*                        zcIsl_Mode = 11 corner cells connect
*                        zcIsl_Mode = 12 corner cells DONT connect
*                        for both modes there is NO connection 
*                        in depth (i.e. EMC<-/->HAC1<-/->HAC2)
*                        since Isles routines recombines this information 
*                        in order to form coneIslands 
*                        from 2-dimentional cellIslands.
*                        IslTyp = 1000 for EMC  cellIslands (2D) 
*                        IslTyp = 2000 for HAC1 cellIslands (2D) 
*                        IslTyp = 3000 for HAC2 cellIslands (2D) 
*
*                        NOTE: it not recomended to use Isles routine with other 
*                              modes than zcIsl_Mode = 11,12 !!!!
*...
      Integer        nIsl
      Integer        NrcIsl,IslTyp
      Real           eIsl,xIsl,yIsl,zIsl,rIsl,emcEIsl,xI1,yI1,zI1
      Common/zIsl02/ nIsl,
     +               NrcIsl( nI_Max),IslTyp( nI_Max),
     +               eIsl(   nI_Max),xIsl(   nI_Max),
     +               yIsl(   nI_Max),zIsl(   nI_Max),
     +               rIsl(   nI_Max),emcEIsl(nI_Max),
     +               xI1(    nI_Max),yI1(    nI_Max),
     +               zI1(    nI_Max)

*..
* 
*..
      Integer        pCellIsl,zIsl_Cell
      Real           ECellIsl,ImbCellIsl
      Common/zIsl03/ pCellIsl(   Nr_cIslMax,nI_Max),
     +               zIsl_Cell(  Nr_cIslMax,nI_Max),
     +               ECellIsl(   Nr_cIslMax,nI_Max),
     +               ImbCellIsl( Nr_cIslMax,nI_Max)

      Integer        I_Caltru
      Common/zIsl04/ I_Caltru(Nr_cMax)

*...
*  (cell/cone)Islands modes....
*  zcIsl_Mode  <==>  1 cells that match at the corner are YES connected
*                      connection in depth (i.e. EMC<->HAC1<->HAC2)
*  zcIsl_Mode  <==>  2 cells that match at the corner are NOT connected
*                      connection in depth (i.e. EMC<->HAC1<->HAC2)
*  zcIsl_Mode  <==> 11 cells that match at the corner are YES connected,
*                      NO connection in depth (2-D mode)
*  zcIsl_Mode  <==> 12 cells that match at the corner are NOT connected, 
*                      NO cennection in depth (2-D mode)
*  zcIsl_Mode  <==> 20 TowerIslands developed by Larry Wai (Columbia)
*
*  NOTE: when using Isles routine (i.e. coneIslands) use zcIsl_Mode=11,12 !!!
*
*  Isles_Debug <==> Debug Flag for more information (lots of output)
*...
      Integer        zcIsl_Mode
      Logical        Isles_Debug
      Common/zIsl05/ Isles_Debug,zcIsl_Mode

