C=======================================================================
C BLOCK CONTROL FOR THE CALCULATION OF THE EDDYVISC
C=======================================================================

      subroutine set_sc1var(u,v,w,sc1,rho,sc1var,sc1dis,
     &                      dip,djp,dkp,dkw,rp,work,offwork)

      implicit none

C PARAMETERS

#        include "mpif.h"
#        include "dim.inc"
#        include "cb.lgridp.inc"
#        include "cb.lgrid.inc"

      real     rho(*)
      real     u(*), v(*), w(*), sc1(*)
      real     sc1var(*), sc1dis(*)
      real     dip(*), djp(*), dkp(*), dkw(*)
      real     rp(*)
      integer  offwork(*)
      real     work(*)

C LOCAL

      integer iblk
      integer pmg_k, pmg_ij
      integer ps_ijk

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

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

        call set_sc1var_block(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &       sc1(ps_ijk),rho(ps_ijk),
     &       dip(pmg_ij),djp(pmg_ij),dkp(pmg_k),dkw(pmg_k),rp(pmg_ij),
     &       sc1var(ps_ijk),sc1dis(ps_ijk))


 10   continue

C-----------------------------------------------------------------------
C PERIODIC BOUNDARIES AND INTERFACES TREAMENTS
C-----------------------------------------------------------------------

      call swap_sc(0,1,sc1var,work,offwork,2)

      return
      end
