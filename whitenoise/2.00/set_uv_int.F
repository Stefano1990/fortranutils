C=======================================================================
C LOOP THROUGH ALL BLOCKS TO SET UP INTERMEDIATE u-VELOCITY
C=======================================================================

      subroutine set_uv_int(iva,u,v,w,p,rho,rhonp1,gam,u_int,hu,
     &  dip,diu,dix,rdix,djp,dju,djv,djx,rdjp,rdju,
     &  didjp,didju,rdidju,
     &  hijp,hiju,hjiu,hkiu,hkju,
     &  rp,ru,dkp,dkw,
     &  tw_sfi,tw_nfi,tw_lfi,tw_rfi,
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

      integer      iva
      real         dip(*), diu(*), dix(*), rdix(*)
      real         djp(*), dju(*), djv(*), djx(*)
      real         rdjp(*), rdju(*)
      real         didjp(*), didju(*), rdidju(*)
      real         hijp(*), hiju(*), hjiu(*), hkiu(*), hkju(*)
      real         rp(*), ru(*)
      
      real         dkp(*), dkw(*)
      
      real         u_int(*)
      real         hu(*)
      real         u(*), v(*), w(*), p(*)
      real         rho(*), gam(*)
      real         rhonp1(*)  
      
      real         tw_sfi(*), tw_nfi(*), tw_lfi(*), tw_rfi(*)
      real         work(*)
      integer      offwork(*)

C LOCAL

      integer      iblk
      integer      nip1, njp1, nkp1
      integer      io, jo, ko, io_ij, jo_ij, ijo(2), ijko(3)
      integer      d1, d2, d12(2,2)
      integer      pm_k, pm_ij, pm_ijk
      integer      ps_ij, ps_ki, ps_ijk
      integer      istart
      logical      walls(6)
      logical      wallwf, wallef, wallsf, wallnf, walllf, wallrf

      data         d12/1,2,2,1/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

C POINTERS FOR MULTIGRID QUANTITIES

        pm_k   = mg_off_1d_k(iblk,1)  + 1
        pm_ij  = mg_off_2d_ij(iblk,1) + 1
        pm_ijk = mg_off(iblk,1)       + 1

C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ij  = sg_off_2d_ij(iblk) + 1
        if (iva .eq. SV_U) then
          ps_ki = sg_off_2d_ki(iblk) + 1
        else
          ps_ki = sg_off_2d_jk(iblk) + 1
        end if
        ps_ijk = sg_off(iblk) + 1

C SET BOUNDARY CONDITION FLAGS

        d1 = d12(1,iva)
        d2 = d12(2,iva)

        if (bl_bc(2*d1-1,iblk) .eq. BCT_WALL) then
          walls(2*d1-1) = .true.
        else
          walls(2*d1-1) = .false.
        end if
        if (bl_bc(2*d1,iblk) .eq. BCT_WALL) then
          walls(2*d1) = .true.
        else
          walls(2*d1) = .false.
        end if
        if (bl_bc(2*d2-1,iblk) .eq. BCT_WALL) then
          walls(2*d2-1) = .true.
        else
          walls(2*d2-1) = .false.
        end if
        if (bl_bc(2*d2,iblk) .eq. BCT_WALL) then
          walls(2*d2) = .true.
        else
          walls(2*d2) = .false.
        end if
        if (bl_bc(5,iblk) .eq. BCT_WALL) then
          walls(5) = .true.
        else
          walls(5) = .false.
        end if
        if (bl_bc(6,iblk) .eq. BCT_WALL) then
          walls(6) = .true.
        else
          walls(6) = .false.
        end if

C SWAP RELEVANT QUANTITIES. NECESSARY WHEN CALLED FOR COMPUTING v_int.

        ijko(1) = 1
        ijko(2) = bl_nijkp1(1,iblk)
        ijko(3) = bl_nijkp1(1,iblk)*bl_nijkp1(2,iblk)

        ijo(1) = 1
        ijo(2) = bl_nijkp1(1,iblk)

        nip1 = bl_nijkp1(d1,iblk)
        njp1 = bl_nijkp1(d2,iblk)
        nkp1 = bl_nijkp1(3 ,iblk)

        io = ijko(d1)
        jo = ijko(d2)
        ko = ijko(3 )

        io_ij = ijo(d1)
        jo_ij = ijo(d2)

        wallwf = walls(2*d1-1)
        wallef = walls(2*d1  )
        wallsf = walls(2*d2-1)
        wallnf = walls(2*d2  )
        walllf = walls(5)
        wallrf = walls(6)

        istart = nijkse(2*d1-1,iblk,1)

        call set_u_int_block(nip1,njp1,nkp1,io,jo,ko,io_ij,jo_ij,dtime,
     &       dip(pm_ij),diu(pm_ij),dix(pm_ij),rdix(pm_ij),djp(pm_ij),
     &       dju(pm_ij),djv(pm_ij),djx(pm_ij),rdjp(pm_ij),rdju(pm_ij),
     &       didjp(pm_ij),didju(pm_ij),rdidju(pm_ij),hijp(pm_ij),
     &       hiju(pm_ij),hjiu(pm_ij),hkiu(pm_ij),hkju(pm_ij),
     &       rp(pm_ij),ru(pm_ij),dkp(pm_k),dkw(pm_k),
     &       u_int(ps_ijk),hu(ps_ijk),u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &       p(pm_ijk),rho(pm_ijk),rhonp1(pm_ijk),gam(ps_ijk),
     &       tw_sfi(ps_ki),tw_nfi(ps_ki),tw_lfi(ps_ij),tw_rfi(ps_ij),
     &       wallwf,wallef,wallsf,wallnf,walllf,wallrf,
     &       istart,iststep,visc,ylaminar,ck,ca,iva)

 10   continue

      call swap_sc(iva,1,u_int,work,2)


      return
      end
