C=======================================================================
C SET SOME FREQUENTLY USED CONSTANTS
C=======================================================================

      subroutine set_const()

      implicit none

#include "dim.inc"

#include "cb.lcontrol.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      cyclic(1,1) = 2
      cyclic(2,1) = 3
      cyclic(1,2) = 3
      cyclic(2,2) = 1
      cyclic(1,3) = 1
      cyclic(2,3) = 2

      cyclic1(1,1) = 2
      cyclic1(2,1) = 3
      cyclic1(1,2) = 1
      cyclic1(2,2) = 3
      cyclic1(1,3) = 1
      cyclic1(2,3) = 2

      return
      end
