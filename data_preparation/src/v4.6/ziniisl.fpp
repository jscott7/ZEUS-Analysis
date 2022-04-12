      Subroutine z_ini_Isles 
      Implicit NONE
*-------------------------------------------------------------------------------
*
*   Called by USER in ZUINIT to initialize parameters for Isles routine
*
*     INPUT:  NONE
*
*     OUTPUT: NONE
*
*   Author:   Gennady Briskin, Tel Aviv University
*   Date:     20-Dec-1996
*
*_______________________________________________________________________________
*
#include "zisles.inc"
*
      Logical First
      DATA    First /.TRUE./
*
      If (First) Then
       First = .FALSE.

       zcIsl_Mode  = 12 
       Isles_Debug = .FALSE.
      EndIf
*
      Write(6,100)Isles_Debug,zcIsl_Mode
*
 100  Format(
     &10x,'*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'/
     &10x,'**                                          **'/
     &10x,'**     Parameters for Isles routine         **'/
     &10x,'**                                          **'/
     &10x,'**       Isles_Debug = ',L8  ,'             **'/
     &10x,'**       zcIsl_Mode  = ',I8  ,'             **'/
     &10x,'**                                          **'/
     &10x,'*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*'/
     &      ) 
      End

