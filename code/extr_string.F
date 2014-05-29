      subroutine extr_string(which_string,string,l_string,string_buffer,
     &  l_buffer,info)
c***********************************************************************
c
c   This subroutine extracts strings from a string buffer
c
c***********************************************************************
c
c   Version   1.0
c   Date      April 19, 1995
c   Computer  JCN PC/486, Sun Sparcstation 10
c
c***********************************************************************

c   Input/Output Parameters

      implicit none
      integer info
      integer which_string, l_string, l_buffer
      character string*(*), string_buffer*(*)

c***********************************************************************
c
c   Variable    i/o   Meaning
c
c***********************************************************************

c   Local Quantities

      integer   l_buf_s, string_seq, string_start, isymbol
      character space, no_more_string, symbol
      logical   is_space, no_more, buf_end, go_on

      data      space/' '/, no_more_string/'!'/

c***********************************************************************
c
c   Explanation
c
c***********************************************************************
c
c   Subroutines Called:
c     none
c
c***********************************************************************

      info = 1

      string_seq   = 0
      string_start = 0
      l_buf_s      = 0
      l_string     = 0

      go_on  = string_seq.ne.which_string.and.l_buf_s.lt.l_buffer

      do while (go_on)

        l_buf_s  = l_buf_s + 1
        symbol   = string_buffer(l_buf_s:l_buf_s)

        isymbol = ichar(symbol)
        is_space = symbol.eq.space.or.isymbol.eq.32.or.isymbol.eq.0
        no_more  = symbol.eq.no_more_string
        buf_end  = l_buf_s.eq.l_buffer
        if(is_space.or.no_more.or.buf_end) then
          if(l_string.ne.0) string_seq = string_seq + 1
          if(string_seq.eq.which_string) then
            go_on = .false.
            info  = 0
          else
            l_string = 0
          end if
        else
          l_string = l_string + 1
          if(l_string.eq.1) string_start = l_buf_s
        end if

        go_on = go_on.and..not.no_more.and.l_buf_s.lt.l_buffer

      end do

      if(buf_end.and..not.no_more) l_string = l_string + 1
c     info != 0 MEANS THAT, FOR INSTANCE, THERE ARE ONLY 5 STRINGS BUT
c     STRINGS AFTER THE 5TH ARE LOOKED FOR
      if(info.eq.0.and.l_string.gt.0) then
        string(1:l_string) =
     &    string_buffer(string_start:string_start+l_string-1)
c     NO STRING FOUND
      else
        info = 2
      end if

c***********************************************************************
c
c   End of subroutine extr_string
c
c***********************************************************************
      end
