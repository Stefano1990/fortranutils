C=======================================================================
C CORRECT THE INTERMEDIATE VELOCITY (ADAMS-BASHFORTH SCHEME)
C=======================================================================

      subroutine uvw_adv(ivar1,diu,rp,dkw,u,u_int,rho,p,wk,offwork,step)

      implicit none

#include "dim.inc"
#include "bc_type.inc"

#include "cb.gridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.comm.inc"
C INPUT/OUTPUT PARAMETERS

      integer ivar1,step
      real    diu(*), rp(*)
      real    dkw(*)
      real    u(*), u_int(*), rho(*), p(*)
      integer offwork(*)
      real    wk(*)

C LOCAL QUANTITIES

      integer iblk
      integer istr, jstr, kstr
      integer nip1, njp1, nkp1
      integer p1, p2, p3, p3_mg
      integer ijkse(6)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

        call set_loop_bound(ivar1,1,1,iblk,iblk,1,ijkse)

        istr = ijkse(1)
        jstr = ijkse(3)
        kstr = ijkse(5)

        p1    = mg_off_1d_k(iblk,1)  + 1
        p2    = mg_off_2d_ij(iblk,1) + 1
        p3_mg = mg_off(iblk,1)       + 1
        p3    = sg_off(iblk)         + 1

        call uvw_adv_block(ivar1,nip1,njp1,nkp1,istr,jstr,kstr,
     &       dtime,diu(p2),rp(p2),dkw(p1),u(p3),u_int(p3),rho(p3_mg),
     &       p(p3_mg),iblkl2g(iblk),step)

 10   continue

      call swap_sc(ivar1,1,u,wk,2)
c     call swap_sc1(ivar1,1,u,wk,offwork,2)

      return
      end
