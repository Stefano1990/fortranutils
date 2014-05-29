C=======================================================================
C MULTIGRID SOLVERS FOR THE SOLUTION OF THE THE PRESSURE EQUATION
C=======================================================================

      subroutine solvers(ap,aw,ae,as,an,al,ar,f,b,residu,ia,ra)

      implicit none

#include "dim.inc"

#include "cb.control.inc"

C INPUT/OUTPUT PARAMETERS

      real     ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*)
      real     f(*), b(*), residu(*)
      real     ra(*)
      integer  ia(*)

C LOCAL QUANTITIES

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      if(ncycle(ivar).eq.0) return

C MULTIGRID

      if(psolver.eq.101) then

        call lmg_solve(ap,aw,ae,as,an,al,ar,f,b,residu,ia,ra)

      else

        write(*,'(''*** Error: No such solver implemented'')')
        stop

      end if

      return
      end
