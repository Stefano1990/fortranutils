C=======================================================================
C RESTRICTION CONTROL
C=======================================================================

      subroutine restriction(ilev,f,f1,work,offwork)

      implicit none

#include "dim.inc"
#include "bc_type.inc"

#include "cb.lgrid.inc"
#include "cb.lcontrol.inc"

      integer ilev
      integer offwork(maxbl,maxmgl)
      real    f(*), f1(*), work(*)

C LOCAL

      integer igridc, igridf
      integer iblk, iff, ifc
      integer nip1c, njp1c, nkp1c, joc, koc
      integer nip1f, njp1f, nkp1f, jof, kof
      integer imult, jmult, kmult

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      igridc = ilev
      igridf = ilev - 1

C-----------------------------------------------------------------------
C LOOP OVER BLOCKS
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        nip1c = mg_nijkp1(1,iblk,igridc)
        njp1c = mg_nijkp1(2,iblk,igridc)
        nkp1c = mg_nijkp1(3,iblk,igridc)
        joc   = nip1c
        koc   = nip1c*njp1c

        nip1f = mg_nijkp1(1,iblk,igridf)
        njp1f = mg_nijkp1(2,iblk,igridf)
        nkp1f = mg_nijkp1(3,iblk,igridf)
        jof   = nip1f
        kof   = nip1f*njp1f

        imult = mg_ijk_mult(1,iblk,igridf)
        jmult = mg_ijk_mult(2,iblk,igridf)
        kmult = mg_ijk_mult(3,iblk,igridf)

        iff = mg_off(iblk,igridf) + 1
        ifc = mg_off(iblk,igridc) + 1

C RESTRICTION

        call rest_block(f(iff),f1(ifc),
     &       nip1c,njp1c,nkp1c,joc,koc,
     &       nip1f,njp1f,nkp1f,jof,kof,
     &       imult,jmult,kmult)

 10   continue

      return
      end
