C=======================================================================
C LOOP THROUGH ALL BLOCKS TO SET VALUES FOR VIRTUAL BOUNDARY POINTS
C=======================================================================

C NOTE: NOT TO BE USED FOR GENERAL BOUNDARY CONDITIONS IN WHICH value IS
C       A FUNCTION OF LOCATION, GRID LEVEL

      subroutine set_virt_bc(iva,ilev,f,fm,rho,gam,
     &  vd_damping,init_vdd,value,
     &  tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &  tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &  tw_wfk,tw_efk,tw_sfk,tw_nfk,
     &  dip,diu,div,dix,djp,dju,djv,djx,ru,rv,dkp,dkw)

      implicit none

#include "dim.inc"
#include "bc_type.inc"
#include "les.inc"
#include "const.inc"
#include "comm.inc"

#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc"

      integer iva, ilev
      real    f(*), fm(*), rho(*), gam(*)
      real    vd_damping(*)
      real    value
      real    tw_sfi(*), tw_nfi(*), tw_lfi(*), tw_rfi(*),
     &        tw_wfj(*), tw_efj(*), tw_lfj(*), tw_rfj(*),
     &        tw_wfk(*), tw_efk(*), tw_sfk(*), tw_nfk(*)
      real    dip(*), diu(*), div(*), dix(*),
     &        djp(*), dju(*), djv(*), djx(*),
     &        ru(*) , rv(*) , dkp(*), dkw(*)
      logical init_vdd(6,*)

C LOCAL

      integer iblk, iside
      integer d1, d2, d3
      integer nip1, njp1, nkp1, nip1_r, njp1_r, nkp1_r
      integer io, jo, ko, kko, jo_ij
      integer k, k1
      integer ijko(3), nijkp1(3)
      integer p1, p2, p3, pij, pki, pkj
      logical walls(6)
      real    rdum(2)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do k = 1,6
        walls(k) = .false.
      end do

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,ilev)
        njp1 = mg_nijkp1(2,iblk,ilev)
        nkp1 = mg_nijkp1(3,iblk,ilev)

        p3 = mg_off(iblk,ilev) + 1

        if (iva.eq.SV_U) then
          do k = 1,6
            init_vdd(k,iblk) = .true.
          end do
        end if

        do iside = p3,p3+nip1*njp1*nkp1-1
          vd_damping(iside) = 1.0
        end do

        do 15 iside = 1,6

          d3 = (iside+1)/2
!wam 12/08/12 zero grad added here for test
          if (bl_bc(iside,iblk).ne.BCT_INTERFACE .and.
     &        bl_bc(iside,iblk).ne.BCT_PERIODIC .and.
!     &        bl_bc(iside,iblk).ne.BCT_ZERO_GRAD .and.
     &        iva.ne.d3) then

            walls(iside) = bl_bc(iside,iblk) .eq. BCT_WALL
            d1 = cyclic1(1,d3)
            d2 = cyclic1(2,d3)

!          write(*,*) iblkl2g(iblk),iside,bl_bc(iside,iblk),walls(iside)

C POINTERS

            p1 = mg_off_1d_k(iblk,ilev) + 1
            p2 = mg_off_2d_ij(iblk,ilev) + 1

            pij = p2
            pki = mg_off_2d_ki(iblk,ilev) + 1
            pkj = mg_off_2d_jk(iblk,ilev) + 1

C ARRANGE LOOPS SO THAT THE STRIDE IS SMALLEST

            ijko(1) = 1
            ijko(2) = nip1
            ijko(3) = nip1*njp1

            nijkp1(1) = nip1
            nijkp1(2) = njp1
            nijkp1(3) = nkp1

            jo_ij = nijkp1(d1)

            io = ijko(d1)
            jo = ijko(d2)
            ko = ijko(d3)

            nip1_r = nijkp1(d1)
            njp1_r = nijkp1(d2)
            nkp1_r = nijkp1(d3)

            if (mod(iside,2) .eq. 1) then
              k   = 1
              k1  = 2
              kko = ko
            else
              k   = nkp1_r
              k1  = nkp1_r
              kko = -ko
            end if

C-----------------------------------------------------------------------
C u-VELOCITY
C-----------------------------------------------------------------------
            if (iva.eq.SV_U .and. ilev.eq.1) then
              if (walls(3)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_sfi(pki),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if(walls(4)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_nfi(pki),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if(walls(5)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_lfi(pij),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if(walls(6)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_rfi(pij),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              end if
C-----------------------------------------------------------------------
C v-VELOCITY
C-----------------------------------------------------------------------
            else if (iva.eq.SV_V .and. ilev.eq.1) then
              if (walls(1)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_wfj(pkj),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if (walls(2)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_efj(pkj),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if (walls(5)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_lfj(pij),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if (walls(6)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_rfj(pij),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              end if
C-----------------------------------------------------------------------
C w-VELOCITY
C-----------------------------------------------------------------------
            else if (iva.eq.SV_W .and. ilev.eq.1) then
              if (walls(1)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_wfk(pkj),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if (walls(2)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_efk(pkj),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if (walls(3)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_sfk(pki),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              else if (walls(4)) then
                call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &               f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &               iside,k,k1,io,jo,ko,kko,
     &               ck,visc,aplus,ylaminar,
     &               walls,init_vdd(1,iblk),tw_nfk(pki),
     &               dip(p2),diu(p2),div(p2),dix(p2),
     &               djp(p2),dju(p2),djv(p2),djx(p2),
     &               ru(p2) ,rv(p2) ,dkp(p1),dkw(p1))
              end if
C-----------------------------------------------------------------------
C OTHER VARIALES, OR VELOCITIES ON NON-WALL BOUNDARY
C-----------------------------------------------------------------------
            else
              call set_virt_bc_block(iva,nip1_r,njp1_r,nkp1_r,
     &             f(p3),fm(p3),rho(p3),gam(p3),vd_damping(p3),value,
     &             iside,k,k1,io,jo,ko,kko,
     &             ck,visc,aplus,ylaminar,
     &             walls,init_vdd(1,iblk),rdum,
     &             rdum,rdum,rdum,rdum,
     &             rdum,rdum,rdum,rdum,
     &             rdum,rdum,rdum,rdum)
            end if

          end if

          walls(iside) = .false.

 15     continue
 10   continue

      return
      end
