!!$ $Id: logging.F90 1759 2013-01-25 11:58:38Z antst $

Module logging
      Implicit None
      Integer, Parameter :: log_stream_number = 1212
      Integer :: have_log_open = 0, log_level = 9, skip_logging = 0
      Character (Len=20) :: logfilename = 'output'
Contains

      Subroutine close_log_file
         Implicit None
         If (have_log_open /= 0) close (log_stream_number)
      End Subroutine close_log_file

      Subroutine flush_log
         Implicit None
         If (have_log_open /= 0) flush (log_stream_number)
      End Subroutine flush_log


      Subroutine do_log (wantlevel, text, pscr)
         Implicit None
         Integer :: wantlevel
         Character (Len=*) :: text
         Logical, Optional :: pscr
!!$
         If (wantlevel > log_level) Return
         If (skip_logging /= 0) Return
         If (have_log_open == 0) Call open_log_file
         Write (log_stream_number, '(A)') text
         flush (log_stream_number)
         If ( .Not. present(pscr)) Return
         If (pscr) Then
            Write (6, '(A)') text
            flush (6)
         End If
      Contains
         Subroutine open_log_file
            Implicit None
            If (have_log_open == 0) Then
               If (skip_logging == 0) Then
                  Open (Unit=log_stream_number, File=logfilename, Action='write')
               Else
                  Open (Unit=log_stream_number, File=logfilename, Action='readwrite')
               End If
               have_log_open = 1
            End If
         End Subroutine open_log_file
      End Subroutine do_log

End Module logging
