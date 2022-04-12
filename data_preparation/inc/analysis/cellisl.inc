

c     Cell Quantities
      integer MaxCell, NCell
      parameter(MaxCell=1000)
      integer
     + CellNr(MaxCell), ! Cell number of ith highest energy Cell
     + CellIsl(MaxCell),! Island number
     + PCIndex(MaxCell),! Index in PC arrays
     + CellList(MaxCell)! List of Cell Numbers Ordered by Island

c     Island Quantities
      integer MaxIsland, nIsland, nPresel(3)
      parameter(MaxIsland=200)
      integer
     + islNCell(MaxIsland),   ! Number of Cells
     + islCellPtr(MaxIsland), ! Offset in CellList for this island
     + islHiCell(MaxIsland),  ! cell number of highest energy cell
     + islHiEMCell(MaxIsland) ! cell number of highest energy EMC cell
      real
     + islEnergy(MaxIsland),  ! energy
     + islHiE(MaxIsland),     ! energy of highest energy cell
     + islFHAC(MaxIsland),    ! HAC fraction (not including HAC0)
     + islFHAC0(MaxIsland),   ! HAC0 fraction
     + islFEMC(MaxIsland),    ! EMC fraction (not including HAC0)
     + islPos(3,MaxIsland),   ! x,y,z position (cm)
     + islSize(3,MaxIsland),  ! x,y,z size (cm)
     + islImbal(-1:1,MaxIsland), ! (EL-ER)/(EL+ER) in each module
     + islHiImbal(MaxIsland), ! imbalance in highest energy cell
     + R21HiMod(MaxIsland),   !
     + fMod(-1:1,MaxIsland),  ! fraction of energy in each module
     + fTow(-1:1,MaxIsland),  ! fraction of energy in each tower
     + fCell(-1:1,MaxIsland)  ! fraction of em energy in each cell
c
      common/CellIsl/ NCell, CellNr, CellIsl, PCIndex, CellList,
     + NIsland, nPresel, islNCell, islCellPtr, islHiCell, islHiEMCell,
     + islEnergy, islHiE, islFHAC, islFHAC0, islFEMC,
     + islPos, islSize, islImbal, islHiImbal,
     + R21HiMod, fMod, fTow, fCell
C
