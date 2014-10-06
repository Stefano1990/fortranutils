C=======================================================================
C LOOP THROUGH ALL BLOCKS TO SET UP INTERMEDIATE w-VELOCITY
C=======================================================================

      subroutine set_w_int(u,v,w,p,rho,rhonp1,gam,w_int,hw,
     &  dip,diu,rdiv,djp,djv,rdju,
     &  didjp,rdidjp,
     &  hkip,hkjp,hkiu,hkjv,
     &  rp,ru,dkp,dkw,
     &  tw_sfk,tw_nfk,tw_wfk,tw_efk,
     &  work,offwork)

      implicit none

C PARAMETERS & INCLUDE FILES

#include "dim.inc"
#include "bc_type.inc"
#include "const.inc"
#include "les.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"

      real         dip(*), diu(*), rdiv(*), djp(*), djv(*), rdju(*)
      real         didjp(*), rdidjp(*)
      real         hkip(*), hkjp(*)
      real         hkiu(*), hkjv(*)
      real         rp(*), ru(*)
      real         dkp(*), dkw(*)
      
      real         w_int(*)
      real         hw(*)
      real         u(*), v(*), w(*), p(*)
      real         rho(*), gam(*)
      real         rhonp1(*) 
      
      real         tw_sfk(*), tw_nfk(*), tw_wfk(*), tw_efk(*)

      real         work(*)
      integer      offwork(*)

C LOCAL

      integer      iblk
      integer      nip1, njp1, nkp1
      integer      io, jo, ko
      integer      pm_k, pm_ij, pm_ijk
      integer      ps_ij, ps_ki, ps_kj, ps_ijk
      integer      kstart
      logical      wallwf, wallef, wallsf, wallnf, walllf, wallrf
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

C SET BOUNDARY CONDITION FLAGS

        if (bl_bc(1,iblk) .eq. BCT_WALL) then
          wallwf = .true.
        else
          wallwf = .false.
        end if
        if (bl_bc(2,iblk) .eq. BCT_WALL) then
          wallef = .true.
        else
          wallef = .false.
        end if
        if (bl_bc(3,iblk) .eq. BCT_WALL) then
          wallsf = .true.
        else
          wallsf = .false.
        end if
        if (bl_bc(4,iblk) .eq. BCT_WALL) then
          wallnf = .true.
        else
          wallnf = .false.
        end if
        if (bl_bc(5,iblk) .eq. BCT_WALL) then
          walllf = .true.
        else
          walllf = .false.
        end if
        if (bl_bc(6,iblk) .eq. BCT_WALL) then
          wallrf = .true.
        else
          wallrf = .false.
        end if

C SWAP RELEVANT QUANTITIES. NECESSARY WHEN CALLED FOR COMPUTING v_int.

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        io = 1
        jo = nip1
        ko = nip1*njp1

        kstart = nijkse(5,iblk,1)

        call set_w_int_block(nip1,njp1,nkp1,io,jo,ko,dtime,
     &       dip(pm_ij),diu(pm_ij),rdiv(pm_ij),djp(pm_ij),djv(pm_ij),
     &       rdju(pm_ij),didjp(pm_ij),rdidjp(pm_ij),hkip(pm_ij),
     &       hkjp(pm_ij),hkiu(pm_ij),hkjv(pm_ij),rp(pm_ij),ru(pm_ij),
     &       dkp(pm_k),dkw(pm_k),
     &       w_int(ps_ijk),hw(ps_ijk),u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &       p(pm_ijk),rho(pm_ijk),rhonp1(pm_ijk),gam(ps_ijk),
     &       tw_sfk(ps_ki),tw_nfk(ps_ki),tw_wfk(ps_kj),tw_efk(ps_kj),
     &       wallwf,wallef,wallsf,wallnf,walllf,wallrf,
     &       kstart,iststep,visc,ylaminar,ck,ca)


 10   continue

      call swap_sc(SV_W,1,w_int,work,2)
c      call swap_sc1(SV_W,1,w_int,work,offwork,2)

C!!!! debug

c     call file_3(20,'w_i',1,w_int,1,2,rdum(1),rdum(2))
c     call file_3(20,'hw',1,hw,1,2,rdum(1),rdum(2))

      return
      end
