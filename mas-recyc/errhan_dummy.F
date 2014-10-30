C=======================================================================
C DUMY ROUTINES FOR SUN'S ERROR HANDLING PROGRAMS (WHEN RUNNING ON OTHER
C MACHINES OTHER THAN SUN, THIS FILE NEEDS TO BE INCLUDED)
C=======================================================================


      integer function ieee_handler()
       ieee_handler=0
      end

      integer function errdivision()
       errdivision=0
      end

      integer function errover()
       errover=0
      end

      integer function errinvalid()
       errinvalid=0
      end

      integer function errunder()
       errunder=0
      end

      integer function errinexact()
       errinexact=0
      end
