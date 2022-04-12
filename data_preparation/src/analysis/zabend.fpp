        SUBROUTINE ZABEND
**********************************************************************
*
*     This routine intercepts ZEBRA ZFATAL calls and attempts to
*   end with some grace, saving histograms, for example.
*
**********************************************************************
      logical    zstop
      save       zstop
      data       zstop / .FALSE. /

      if (zstop) then
C        Protect against reentrancy
C        --------------------------
         STOP
      else
         zstop = .TRUE.
         WRITE (6, '(''1  Job Terminated with a ZEBRA Error'')' )

         CALL ZUTERM
      endif

      RETURN
      END









