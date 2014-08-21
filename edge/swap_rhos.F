C=======================================================================
C BLOCK CONTROL FOR THE CALCULATION OF THE dp(n+1)/dt
C=======================================================================

      subroutine swap_rhos(rho,rhonp1,rhonm1)

      implicit none

C PARAMETERS

#include "dim.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"

      real     rho(*), rhonp1(*), rhonm1(*)

C LOCAL

      integer      iblk
      integer      nip1, njp1, nkp1
      integer p1, p2, p3
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C CALCULATE EDDY VISCOSITY BLOCK BY BLOCK
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

C POINTERS FOR MULTIGRID RELATED QUANTITIES

        p1  = mg_off_1d_k(iblk,1)  + 1
        p2  = mg_off_2d_ij(iblk,1) + 1
        p3  = mg_off(iblk,1)       + 1

        call swap_rhos_block(rho(p3),rhonp1(p3),rhonm1(p3),
     &                       nip1, njp1, nkp1)

 10   continue


      return
      end
