C=======================================================================
C SOLVER CHOICE FOR THE SOLUTION ON THE COARSEST GRID
C 1. - CURRENTLY ONLY THE SMOOTHER IS USED AS THE SOLVER
C=======================================================================


      subroutine solversc(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,niters,
     &  omega,soltol,work,offwork,resnrm,rhsnrm,errnrm,varnrm)

      implicit none

#include "dim.inc"

#include "cb.control.inc"

C INPUT/OUTPUT PARAMETERS

      integer ilev, niters
      real    ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*)
      real    f(*), b(*), residu(*)
      real    omega, soltol
      real    resnrm, rhsnrm, errnrm, varnrm
      real    work(*)
      integer offwork(*)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C THAT meth_csolve=0 MEANS USING THE SMOOTHER AS THE SOLVER ON THE
C COARSEST GRID

      if(meth_csolve.eq.0) then

        call smooth(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,niters,
     &       omega,tolc,work,offwork,resnrm,rhsnrm,errnrm,varnrm)

      else

C OTHER SOLVERS ARE NOT YET IMPLEMENTED

        write(*,'(''*** Error: '',
     &            ''No such solution method on the coarsest grid'')')
        stop

      end if

      return
      end
