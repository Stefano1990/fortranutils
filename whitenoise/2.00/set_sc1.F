C=======================================================================
C LOOP THROUGH ALL BLOCKS TO SET UP SCALAR
C=======================================================================

      subroutine set_sc1(u,v,w,p,sc1,rho,gamsc1,gamlam,hsc1,
     &  dip,diu,rdiv,djp,djv,rdju,
     &  didjp,rdidjp,
     &  hkip,hkjp,hkiu,hkjv,
     &  rp,dkp,dkw,
     &  work,offwork)

      implicit none

C PARAMETERS & INCLUDE FILES

#include "mpif.h"
#include "dim.inc"
#include "bc_type.inc"
#include "const.inc"
#include "les.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"

      real         dip(*), diu(*), rdiv(*), djp(*), djv(*), rdju(*)
      real         didjp(*), rdidjp(*)
      real         hkip(*), hkjp(*)
      real         hkiu(*), hkjv(*)
      real         rp(*)
      real         dkp(*), dkw(*)
      
      real         hsc1(*)
      real         u(*), v(*), w(*), p(*), sc1(*)
      real         rho(*), gamsc1(*),gamlam(*)
      
      real         work(*)
      integer      offwork(*)

C LOCAL

      integer      iblk
      integer      nip1, njp1, nkp1
      integer      io, jo, ko
      integer      pm_k, pm_ij, pm_ijk
      integer      ps_ij, ps_ki, ps_kj, ps_ijk
      integer      kstart
      real         rdum(2)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

C POINTERS FOR MULTIGRID QUANTITIES

        pm_k   = mg_off_1d_k(iblk,1)  + 1
        pm_ij  = mg_off_2d_ij(iblk,1) + 1
        pm_ijk = mg_off(iblk,1)       + 1

C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ij  = sg_off_2d_ij(iblk) + 1
        ps_ki  = sg_off_2d_ki(iblk) + 1
        ps_kj  = sg_off_2d_jk(iblk) + 1
        ps_ijk = sg_off(iblk)       + 1

C SWAP RELEVANT QUANTITIES. NECESSARY WHEN CALLED FOR COMPUTING v_int.

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        io = 1
        jo = nip1
        ko = nip1*njp1

        kstart = nijkse(5,iblk,1)

        call set_sc1_block(nip1,njp1,nkp1,io,jo,ko,dtime,
     &       dip(pm_ij),diu(pm_ij),rdiv(pm_ij),djp(pm_ij),djv(pm_ij),
     &       rdju(pm_ij),didjp(pm_ij),rdidjp(pm_ij),rp(pm_ij),
     &       dkp(pm_k),dkw(pm_k),
     &       u(ps_ijk),v(ps_ijk),w(ps_ijk),p(pm_ijk),sc1(ps_ijk),
     &       hsc1(ps_ijk),rho(pm_ijk),gamsc1(ps_ijk),gamlam(ps_ijk),
     &       visc,kstart,iststep,iblkl2g(iblk))

 10   continue

      call swap_sc(SV_SC1,1,sc1,work,2)
c      call swap_sc1(SV_SC1,1,sc1,work,offwork,2)

      return
      end
