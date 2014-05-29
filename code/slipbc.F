!------------------------------------
!test!
      subroutine slipbc(x,y,z,dip,diu,djp,djv,dkp,
     &             u,v,v_int,w,sc1,rho,p,step)

       implicit none
! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
#include "data_type.inc"

      integer step,iblk
      real  u(*), v(*), w(*), sc1(*), rho(*), p(*), v_int(*)
      real  dip(*), diu(*), djp(*),
     &      djv(*)
      real  dkp(*), x(*),y(*), z(*)

      integer  pmg_k, pmg_ij,p_xy
      integer  ps_ijk, ps_wijk, ps_wij


!        write(*,*) ' in slipbc'

      do iblk = 1,numbl

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1
        p_xy = mg_off_2d_xy(iblk,1) + 1
C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1

       ps_wijk = 1 + ((iblk-1)*nip1*njp1*nkp1)
       ps_wij = 1 + (iblk-1)*nip1*njp1

      call slipbc_block(u(ps_ijk),v(ps_ijk),v_int(ps_ijk),w(ps_ijk),
     &      rho(ps_ijk),sc1(ps_ijk), p(ps_ijk),dip(pmg_ij),
     &      diu(pmg_ij),djp(pmg_ij),
     &      djv(pmg_ij),dkp(pmg_k),nip1,njp1,nkp1,iblk,step)



      enddo

      return
      end

!------------------------------------------------


      subroutine slipbc_block(u,v,vint,w,rho,
     &      sc1,p,dip,diu,djp,
     &      djv,dkp,nip1,njp1,nkp1,ibl,step)

      implicit none
C PARAMETERS

#include "dim.inc"
#include "bc_type.inc"
#include "const.inc"
#include "les.inc"
#include "comm.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"

      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     rho(nip1*njp1*nkp1), sc1(nip1*njp1*nkp1) 
      real     p(nip1*njp1*nkp1),vint(nip1*njp1*nkp1)
      real     dip(nip1,njp1), diu(nip1,njp1), djp(nip1,njp1),
     &         djv(nip1,njp1)
      real     dkp(nkp1), rhou, ddtime
!local
      integer  i, j, k, jk, ijk, ijkw, ijke, ijks, ijkn, ijkl
      integer  io,jo,ko,ibl,nip1,njp1,nkp1,step


      ddtime = 1.5*dtime
               
               jo = nip1
               ko = nip1*njp1
               if (bl_bc(3,ibl) .eq. BCT_TRACTION) then


!               write(*,*)'traction bc on side 3', iblkl2g(ibl),imax,jmax
                do i = 1, nip1 
                   do k = 2, nkp1-1
                      j=2
                      ijk = i + jo*(j-1) + ko*(k-1) 
                      ijks = ijk - jo
                      ijkn = ijk + jo
                      ijkw = ijk-1
                      ijkl = ijk - ko
                      sc1(ijks) = sc1(ijk)
                      rho(ijks) = rho(ijk)
                      rhou = 0.5*(rho(ijk) + rho(ijks))
                      u(ijks) = u(ijk)
!                      v(ijk) = 0.0
!                      v(ijk) = 0.0 - ddtime*((2.0*p(ijk) 
!     &                         - 2.0*0.0)/djp(i,1)) / rhou
                      v(ijk) = vint(ijk) - ddtime*((2.0*p(ijk)
     &                         - 2.0*0.0)/djp(i,1)) / rhou
                      v(ijks) = v(ijk)
                      w(ijks) = w(ijk)
                   end do
                end do

             end if

             if (bl_bc(4,ibl) .eq. BCT_TRACTION) then

!               write(*,*) 'traction bc on side 4', iblkl2g(ibl)
                do i = 1, nip1
                   do k = 2, nkp1-1
                      j = njp1
                      ijk = i + jo*(j-2) + ko*(k-1)
                      ijkn = ijk + jo 
                      ijkw = ijk-1
                      ijks = ijk - jo
                      ijkl = ijk - ko
                      sc1(ijkn) = sc1(ijk)
                      rho(ijkn) = rho(ijk)
                      rhou = 0.5*(rho(ijk) + rho(ijkn))
                      u(ijkn) = u(ijk)
                      v(ijkn) = vint(ijkn) - ddtime*((-2.0*p(ijk) 
     &                         + 2.0*0.0)/djp(i,njp1-1)) / rhou
!                      v(ijkn) = 0.0 - ddtime*((-2.0*p(ijk)
!     &                         + 2.0*0.0)/djp(i,njp1-1)) / rhou
                      w(ijkn) = w(ijk)
                   end do
                end do

             end if


               if (bl_bc(3,ibl) .eq. BCT_ZERO_GRAD) then


!               write(*,*)'traction bc on side 3', iblkl2g(ibl),imax,jmax
                do i = 1, nip1 
                   do k = 2, nkp1-1
                      j=2
                      ijk = i + jo*(j-1) + ko*(k-1) 
                      ijks = ijk - jo
                      ijkn = ijk + jo
                      ijkw = ijk-1
                      ijkl = ijk - ko
                      sc1(ijks) = sc1(ijk)
                      rho(ijks) = rho(ijk)
                      u(ijks) = u(ijk)
                      v(ijk) = 0.0
                      v(ijks) = v(ijk)
                      w(ijks) = w(ijk)
                   end do
                end do

             end if

             if (bl_bc(4,ibl) .eq. BCT_ZERO_GRAD) then

!               write(*,*) 'traction bc on side 4', iblkl2g(ibl)
                do i = 1, nip1
                   do k = 2, nkp1-1
                      j = njp1
                      ijk = i + jo*(j-2) + ko*(k-1)
                      ijkn = ijk + jo 
                      ijkw = ijk-1
                      ijks = ijk - jo
                      ijkl = ijk - ko
                      sc1(ijkn) = sc1(ijk)
                      rho(ijkn) = rho(ijk)
                      u(ijkn) = u(ijk)
                      v(ijkn) = 0.0
                      w(ijkn) = w(ijk)
                   end do
                end do

             end if




       return

      end

