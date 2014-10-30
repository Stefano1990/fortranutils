C=======================================================================
C EXTRAPOLATE PRESSURE TO THE BOUNDARY POINTS
C=======================================================================

      subroutine extrap_p(p,work,offwork)

      implicit none

#include "dim.inc"
#include "les.inc"

#include "cb.lgrid.inc"

C INPUT/OUTPUT PARAMETERS

      real    p(*)
      integer offwork(*)
      real    work(*)

C LOCAL QUANTITIES

      integer iblk, p3
      integer nip1, njp1, nkp1

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C CALL TO extrap_p_block TO EXTRAPOLATE THE PRESSURE IN A BLOCK

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

        p3 = mg_off(iblk,1) + 1

        call extrap_p_block(p(p3),nip1,njp1,nkp1,bl_bc(1,iblk))

 10   continue

C SWAP INFORMATION FOR INTERFACES AND PERIODODIC BOUNDARIES

      call swap_sc(SV_P,1,p,work,1)
c     call swap_sc1(SV_P,1,p,work,offwork,1)

      return
      end
