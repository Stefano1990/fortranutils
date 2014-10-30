C=======================================================================
C BLOCK CONTROL FOR THE CALCULATION OF
C=======================================================================

      subroutine mass_change(rho,rhonp1,mass_chg,rdidjp,dkw)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS

#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

C COMMON BLOCKS
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

C INPUT/OUTPUT PARAMETERS

      real     rho(*), rhonp1(*)
      real     mass_chg, mass_chg_blk

      real     rdidjp(*) 
      real     dkw(*)   

C LOCAL

      integer iblk
      integer pmg_k, pmg_ij
      integer ps_ijk

	real red

	integer ierr

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------
       mass_chg = 0.0

C-----------------------------------------------------------------------
C CALCULATE EDDY VISCOSITY BLOCK BY BLOCK
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

C POINTERS FOR MULTIGRID RELATED QUANTITIES

        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1

C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1

        call mass_change_block(rho(ps_ijk),rhonp1(ps_ijk),
     &  mass_chg_blk,rdidjp(pmg_ij),dkw(pmg_k),dtime)

        mass_chg = mass_chg + mass_chg_blk

 10   continue



        red = mass_chg

        call MPI_ALLREDUCE(red, mass_chg, 1, REAL_TYPE_MPI,

     &  MPI_SUM, MPI_COMM_WORLD, ierr)


      return
      end
