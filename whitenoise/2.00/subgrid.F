!-------------------------------------------------------
! Controlling routine for calculation of eddy viscosity
!-------------------------------------------------------

      subroutine subgrid(gam,gamsc1,u,v,w,sc1,rho,sc1var,sc1dis,
     &  vd_damping,dip,diu,djp,djv,dkp,dkw,
     &  hijp,hjip,hkip,hkjp,rp,work,offwork)

      implicit none

C PARAMETERS

#include "dim.inc"

#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"

      real     gam(*), rho(*), vd_damping(*), gamsc1(*)
      real     u(*), v(*), w(*), sc1(*)
      real     sc1var(*), sc1dis(*)
      real     dip(*), diu(*), djp(*), djv(*), dkp(*), dkw(*)
      real     hijp(*), hjip(*), hkip(*), hkjp(*)
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

      do iblk = 1,numbl

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

C POINTERS FOR MULTIGRID RELATED QUANTITIES

        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1

C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1

        if(model .eq. Smagorinsky) then

        call subgrid_smagorinsky(gam(ps_ijk),gamsc1(ps_ijk),u(ps_ijk),
     &       v(ps_ijk),w(ps_ijk),sc1(ps_ijk),rho(ps_ijk),
     &       vd_damping(ps_ijk),
     &       dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),djv(pmg_ij),dkp(pmg_k),
     &       dkw(pmg_k),
     &       hijp(pmg_ij),hjip(pmg_ij),hkip(pmg_ij),hkjp(pmg_ij),
     &       rp(pmg_ij),bl_bc(1,iblk),sc1var(ps_ijk),
     &       sc1dis(ps_ijk))
     
         endif
         
       if(model .eq. Wale) then
         call subgrid_wale(gam(ps_ijk),gamsc1(ps_ijk),u(ps_ijk),
     &       v(ps_ijk),w(ps_ijk),sc1(ps_ijk),rho(ps_ijk),
     &       vd_damping(ps_ijk),
     &       dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),djv(pmg_ij),dkp(pmg_k),
     &       dkw(pmg_k),
     &       hijp(pmg_ij),hjip(pmg_ij),hkip(pmg_ij),hkjp(pmg_ij),
     &       rp(pmg_ij),bl_bc(1,iblk),sc1var(ps_ijk),
     &       sc1dis(ps_ijk))        
         
        endif
        
        if(model .eq. None) then
        call subgrid_nosgs(gam(ps_ijk),gamsc1(ps_ijk),u(ps_ijk),
     &       v(ps_ijk),w(ps_ijk),sc1(ps_ijk),rho(ps_ijk),
     &       vd_damping(ps_ijk),
     &       dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),djv(pmg_ij),dkp(pmg_k),
     &       dkw(pmg_k),
     &       hijp(pmg_ij),hjip(pmg_ij),hkip(pmg_ij),hkjp(pmg_ij),
     &       rp(pmg_ij),bl_bc(1,iblk),sc1var(ps_ijk),
     &       sc1dis(ps_ijk))     
        endif
     
      enddo

C-----------------------------------------------------------------------
C PERIODIC BOUNDARIES AND INTERFACES TREAMENTS
C-----------------------------------------------------------------------

      call swap_sc(0,1,gam,work,2)
c      call swap_sc1(0,1,gam,work,offwork,2)

      call swap_sc(0,1,gamsc1,work,2)
c      call swap_sc1(0,1,gamsc1,work,offwork,2)

      return
      end
