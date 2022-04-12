*Start_Documentation
*********************************************************************
*     Common Blocks for "toolbox" (PATCH TOOLNPA, TOOLBOX)
*********************************************************************
*End_Documentation
      COMMON /LEPKIN/ ELBEAM(4),PRBEAM(4),ELSCAP(4),QVGAPP(4)
     +         ,BREMSG(4),QVGTRU(4)
     +         ,XBJAPP,XBJTRU,YPSAPP,YPSTRU,QSQAPP,QSQTRU
     +         ,IDELBM,IDELSC,IDBREM,ICHNN
      REAL*8  ELBEAM,PRBEAM,ELSCAP,QVGAPP,BREMSG,QVGTRU
     +         ,XBJAPP,XBJTRU,YPSAPP,YPSTRU,QSQAPP,QSQTRU
      INTEGER  IDELBM,IDELSC,IDBREM,ICHNN
*
