        Subroutine XBacSum(Vertex,Imod,Ebac,Etbac,pxbac,pybac,pzbac)
*       ************************************************************
*
C *****************************************************************************
C#
C#        Subroutine XBacSum(Vertex,Imod,Ebac,Etbac,pxbac,pybac,pzbac)
C#        ============================================================
C#
C@#   XBacSum:       Calculate total BAC energy and momentum
C#
C#  Purpose        : Calculate total BAC energy, transvers energy and
C#                   momentum vector for energy clusters found in BAC.
C#                   All condensates (Imod=1) or only those matched to CAL
C#                   jets (Imod=2) are used. BAC deposits are corrected
C#                   for energy losses between CAL and BAC.
C#
C#  This version   : In new BAC reconstruction clusters can be created
C#                   also based on wire or strip information (without pad
C#                   energy). The routine however considers only clusters
C#                   containing pads ( XJetEt_ClClass = 0 or 1 ).
C#
C#  Arguments      : Vertex(3) - vertex position as used for jet calculation
C#                   Imod      - flag use to choose operation mode:
C#                               = 1 use all BAC condensates
C#                               = 2 use only BAC condensates mached to
C#                                   CAL jets (i.e. linked in XMatEt table)
C#
C#  Output value   : Ebac      - total energy of BAC condensates
C#                   Etbac     - total transverse energy
C#                   pxbac,pybac,pzbac
C#                             - total momentum vector in BAC
C#              !!!  All values include correction for energy losses
C#                   in dead material between CAL and BAC
C#
C# Implicit inputs : XJetEt and XMipEt (BAC condensate tables)
C#                   XMatEt - table used for BAC to CAL matching
C#
C# Implicit outputs: None
C#
C# Side effects    : None
C#
C# Examples        : see BAC WWW page or
C#                   ~zarnecki/bac/data95/export/export.readme
C#
C# Author          : Aleksander Zarnecki, February 1996.
C#                   zarnecki@vxdesy
C#
C *****************************************************************************
C     Date   ; Name            ; Description
C -----------+-----------------+-----------------------------------------------
C 28.11.1996 |A.F.Zarnecki     | Corrected for new DDL definition of XJetEt.
C            |                 | Only pad clusters are concidered.
C -----------+-----------------+-----------------------------------------------
C *****************************************************************************
*
        Implicit NONE
*
      Real    Vertex(3)
      Integer Imod
      Real    Ebac,Etbac,pxbac,pybac,pzbac
      LOGICAL debug
      DATA    debug /.false./
*
#include "partap.inc"
#include "xjetet.inc"
#include "xmipet.inc"
#include "xmatet.inc"
*
       Integer Iflag
       Integer Ncluj,Nclum,Nclu,Nmat,Iclu
       Real Rxy,PhiBac,ThetaBac
       Real Cbac,cdef,Eleak
       Parameter (Cdef=1.40)
*
       Logical First,RelFlag
       Data First /.TRUE./
*
*   Test printout flag:
*
       Logical Lprint
       Integer Iprint,Nprint
       Data Iprint /0/
        Parameter (Nprint=0)

      IF (debug) WRITE(*,*) 'start XBACSUM'
*
*-----------------------------------------------------------------------
*
*  For first event
*
       IF (First) Then
*
         Write(6,100)
 100     Format(5X,50('*')/
     +       5X,'*',48X,'*'/
     +       5X,'*    XBacSum - sum jet energy leakages found     *'/
     +       5X,'*             in BAC calorimeter                 *'/
     +       5X,'*',48X,'*'/
     +       5X,'*            Version 2 - November 1996           *'/
     +       5X,'*',48X,'*'/
     +       5X,50('*'))

          First=.false.
*
          EndIf
*
*-----------------------------------------------------------------------
*
*  Check mode flag
*
       Iflag=Imod
*
       If(Iflag.NE.1  .AND. Iflag.NE.2)Then
          Write(6,101)Iflag
  101     Format(5X,'XBacSum: ERROR ! Wrong mode Imod = ',I10/
     +       14X,'Default mode Imod = 1 used (all BAC clusters)')
          Iflag=1
          EndIf
*
C----------------------------------------------------------------
C
C  Test printout flag
C
       Lprint = (Iprint.LT.Nprint)
C
*-----------------------------------------------------------------------
*
       Ebac  = 0.
       Etbac = 0.
       pxbac = 0.
       pybac = 0.
       pzbac = 0.
*
*-----------------------------------------------------------------------
*
*  Check number of BAC clusters, return if 0.
*
       Ncluj = CouTab(XJetEt)
       Nclum = CouTab(XMipEt)
       Nclu  = Ncluj+Nclum
*
       If(Nclu.LE.0)Return
*
*-----------------------------------------------------------------------
*
*  Check XMatEt table (should not be empty, if Imod=2)
*
       Nmat = CouTab(XMatEt)

       If(Iflag.EQ.2  .AND. Nmat.LE.0)Then
          Write(6,102)
  102     Format(5X,'XBacSum: Warrning ! XMatEt table empty '/
     +       14X,'Mode Imod = 1 used (all BAC clusters)')
          Iflag=1
          EndIf
*
*-----------------------------------------------------------------------
*
*  Loop over BAC clusters
*
*   BAC "jets":
*
       Do Iclu = 1,Ncluj

         Call FetTab(XJetEt,ID,Iclu)
*
* Take only pad clusters (or output of old reconstruction - 'unknown')
*
         If(XJetEt_ClClass.LE.1 .OR. XJetEt_ClClass.EQ.INULL)Then

            XMatEt_EUCAL = 0.
            XMatEt_EBAC  = XjetEt_XEnDep
            XMatEt_ESum  = Cdef*XjetEt_XEnDep
*
*  Check matching with CAL
*
C ---       AQ hack to avoid error statements in log file ---
            IF (XJetEt_XMatEt.GT.0.AND.
     &          XJetEt_XMatEt.LE.COUTAB(XMatEt)) THEN
              Call NaTRel(XJetEt,XJetEt_XMatEt,XMatEt,RelFlag)
            ELSE
              Iflag = 0
            ENDIF

            If(XMatEt_EUCAL.GT.0. .OR. Iflag.EQ.1)Then

               PhiBac=atan2(XJetEt_XYZCOG(2)-Vertex(2),
     +                              XJetEt_XYZCOG(1)-Vertex(1))
               Rxy=sqrt((XJetEt_XYZCOG(1)-Vertex(1))**2+
     +                             (XJetEt_XYZCOG(2)-Vertex(2))**2)
               ThetaBac=atan2(Rxy,XJetEt_XYZCOG(3)-Vertex(3))

               Eleak = XMatEt_ESum - XMatEt_EUCAL
               Cbac  = Eleak/XMatEt_EBAC

               Eleak = Cbac*XjetEt_XEnDep

               Ebac  = Ebac + Eleak
               pzbac = pzbac + Eleak * cos(Thetabac)

               Eleak = Eleak * sin(ThetaBac)

               Etbac = Etbac + Eleak
               pxbac = pxbac + Eleak * cos(Phibac)
               pybac = pybac + Eleak * sin(Phibac)

               EndIf

            EndIf

         Enddo
*
*    BAC "mips" - from old reconstruction
*
       Do Iclu=1,Nclum

         Call FetTab(XMipEt,ID,Iclu)
*
*  Check if this cluster is not ambigious,
*  if so: skip it
*
         Call NaTRel(XMipEt,XMipEt_XJetEt,XJetEt,RelFlag)

         If( .NOT.RelFlag ) Then
*
            XMatEt_EUCAL=0.
            XMatEt_EBAC = XMipEt_XEnDep
            XMatEt_ESum = Cbac*XMipEt_XEnDep
*
*  Check matching with CAL
*
            Call NaTRel(XMipEt,XMipEt_XMatEt,XMatEt,RelFlag)

            If (XMatEt_EUCAL.GT.0. .OR. Iflag.EQ.1) Then

               PhiBac=atan2(Xmipet_XYZFPO(2)-Vertex(2),
     +                              Xmipet_XYZFPO(1)-Vertex(1))
               Rxy=sqrt((Xmipet_XYZFPO(1)-Vertex(1))**2+
     +                             (Xmipet_XYZFPO(2)-Vertex(2))**2)
               ThetaBac=atan2(Rxy,Xmipet_XYZFPO(3)-Vertex(3))

               Eleak = XMatEt_ESum - XMatEt_EUCAL
               Cbac  = Eleak/XMatEt_EBAC

               Eleak = Cbac*XMipEt_XEnDep

               Ebac  = Ebac + Eleak
               pzbac = pzbac + Eleak * cos(Thetabac)

               Eleak = Eleak * sin(ThetaBac)

               Etbac = Etbac + Eleak
               pxbac = pxbac + Eleak * cos(Phibac)
               pybac = pybac + Eleak * sin(Phibac)

               EndIf

            EndIf

         Enddo
*
*-----------------------------------------------------------------------
*
*  Test printout
*
       If(Lprint)Then

          If(Iflag.EQ.1)Then
             Write(6,*)'    XBacSum sum over all BAC condensates:'
          Else
             Write(6,*)'    XBacSum sum over CAL - BAC matches:'
             EndIf

          Write(6,103)Ebac,Etbac,pxbac,pybac,pzbac

 103      Format(5X,'E = ',F6.2,', Et = ',F6.2,', px = ',F6.2,
     +      ', py = ',F6.2,', pz = ',F6.2/)

        Iprint=Iprint+1

          EndIf


      IF (debug) WRITE(*,*) 'finish XBACSUM'
*
*-----------------------------------------------------------------------
*
       Return
       End
