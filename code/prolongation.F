C=======================================================================
C PROLONGATION CONTROL FOR MULTIBLOCK
C=======================================================================

      subroutine prolongation(ilev,f,work,offwork)

      implicit none

#include "mpif.h"
#include "dim.inc"
#include "data_type.inc"
#include "comm.inc"
#include "les.inc"

#include "cb.lgrid.inc"
#include "cb.comm.inc"

C INPUT/OUTPUT

      integer ilev
      real    f(*), work(*)
      integer offwork(*)

C LOCAL

      integer iblk, iblk0
      integer igridc, igridf
      integer nip1c, njp1c, nkp1c, nip1f, njp1f, nkp1f
      integer imult, jmult, kmult
      integer iff, ifc

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C NUMBERS OF ACTIVE PES ON THE COARSE AND FINE GRIDS

C - l STARTS FROM 1 FOR THE COARSEST GRID. THE INPUT VALUE OF IT HERE
C   IS FOR THE FINER GRID OF THE TWO GRID LEVELS INVOLVED IN PROLONGA-
C   TION. igrid STARTS FROM 1 FOR THE FINEST GRID LEVEL

      igridf = ilev
      igridc = igridf + 1

      do 10 iblk = 1,numbl

        nip1c = mg_nijkp1(1,iblk,igridc)
        njp1c = mg_nijkp1(2,iblk,igridc)
        nkp1c = mg_nijkp1(3,iblk,igridc)
        nip1f = mg_nijkp1(1,iblk,igridf)
        njp1f = mg_nijkp1(2,iblk,igridf)
        nkp1f = mg_nijkp1(3,iblk,igridf)

        imult = mg_ijk_mult(1,iblk,igridf)
        jmult = mg_ijk_mult(2,iblk,igridf)
        kmult = mg_ijk_mult(3,iblk,igridf)

C-----------------------------------------------------------------------
C PROLONGATION
C-----------------------------------------------------------------------

        iff = mg_off(iblk,igridf) + 1
        ifc = mg_off(iblk,igridc) + 1

c       iblk0 = iblkl2g(iblk)
        call prol_block(f(iff),f(ifc),work,
     &       nip1c,njp1c,nkp1c,nip1f,njp1f,nkp1f,
     &       imult,jmult,kmult,bl_bc(1,iblk))

 10   continue

      call swap_sc(SV_P,igridf,f,work,1)
c     call swap_sc1(SV_P,igridf,f,work,offwork,1)

      return
      end
