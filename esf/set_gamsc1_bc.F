C=======================================================================
C LOOP THROUGH ALL BLOCKS TO SET UP BOUNDARY VALUES FOR THE SCALAR
C=======================================================================

      subroutine set_gamsc1_bc(sc1,gamsc1,gamlam,diu,djv,dkp,rp)

      implicit none

C PARAMETERS & INCLUDE FILES

#include "dim.inc"
#include "bc_type.inc"
#include "const.inc"
#include "les.inc"
#include "comm.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"
#include "tvdgrad.inc"

      real         sc1(*), gamsc1(*), gamlam(*)
      real         diu(*), djv(*), dkp(*), rp(*)
      real         gamtmp
      
C LOCAL

      integer      iblk, iside, i, j, k, ij, ik, jk, jk1
      integer      ij1
      integer      nip1, njp1, nkp1, ni, nj, nk
      integer      io, jo, ko
      integer      ps_ijk, pm_ij, pm_k, ps_ijk_last
      integer      ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr, ijk2, ijk1
      integer      ijk3

      real         sclam
C-----------------------------------------------------------------------

      sclam = 0.72

      do iblk = 1,numbl

        ps_ijk = sg_off(iblk)

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        jo = nip1
        ko = nip1*njp1

        do i = 1, nip1
           do j = 1, njp1
              do k = 1, nkp1
                 ijk = i + jo*(j-1) + ko*(k-1) + ps_ijk
                 gamlam(ijk) = visc / sclam
              enddo
           enddo
        enddo

      enddo

      do 10 iblk = 1,numbl

        ps_ijk = sg_off(iblk)
        pm_k = mg_off_1d_k(iblk,1)
        pm_ij = mg_off_2d_ij(iblk,1)

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        ni = nip1 - 1
        nj = njp1 - 1
        nk = nkp1 - 1

        io = 1
        jo = nip1
        ko = nip1*njp1

        do j = 2, nj
           do k = 2, nk
              jk = jo*(j-1) + ko*(k-1) + ps_ijk 
              iside = 1
              ij1 = 3 + jo*(j-1) + pm_ij
              ijkw = 1 + jk
              ijk = 2 + jk
              ijk1 = ijk + 1
              ijk2 = ni + jk
              ijk3 = nip1 + jk
c              dfdxw(iblk,j,k) = (sc1(ijk1) - sc1(ijk)) / diu(3,j)
              dfdxw(iblk,j,k) = (sc1(ijk1) - sc1(ijk)) / diu(ij1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 gamsc1(ijkw) = gamsc1(ijk2) 
              else if (bl_bc(iside,iblk) .eq. BCT_WALL) then
                 gamsc1(ijkw) = -gamsc1(ijk)
                 gamlam(ijkw) = -gamlam(ijk)
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE) then
                 gamsc1(ijkw) = gamsc1(ijk)
              endif
              iside = 2
              ij1 = ni + jo*(j-1) + pm_ij
              ijke = nip1 + jk
              ijk = ijke - 1
              ijk1 = ijk - 1
              ijk2 = 2 + jk
c              dfdxe(iblk,j,k) = (sc1(ijk) - sc1(ijk1)) / diu(ni,j)
              dfdxe(iblk,j,k) = (sc1(ijk) - sc1(ijk1)) / diu(ij1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 gamsc1(ijke) = gamsc1(ijk2) 
              else if (bl_bc(iside,iblk) .eq. BCT_WALL) then
                 gamsc1(ijke) = -gamsc1(ijk)
                 gamlam(ijke) = -gamlam(ijk)
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE) then
                 gamsc1(ijke) = gamsc1(ijk)
              endif
           enddo
        enddo

        do i = 2, ni
           do k = 2, nk
              ik = i + ko*(k-1) + ps_ijk 
              iside = 3
              ij1 = i + jo*(3-1) + pm_ij
              ijks = ik + jo*(1-1)
              ijk = ik + jo*(2-1)
              ijk1 = ijk + jo
              ijk2 = ik + jo*(nj-1)
c              dfdys(iblk,i,k) = (sc1(ijk1) - sc1(ijk)) / djv(i,3)
              dfdys(iblk,i,k) = (sc1(ijk1) - sc1(ijk)) / djv(ij1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 gamsc1(ijks) = gamsc1(ijk2) 
              else if (bl_bc(iside,iblk) .eq. BCT_WALL) then
                 gamsc1(ijks) = -gamsc1(ijk)
                 gamlam(ijks) = -gamlam(ijk)
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE) then
                 gamsc1(ijks) = gamsc1(ijk)
              endif
              iside = 4
              ij1 = i + jo*(nj-1) + pm_ij
              ijkn = ik + jo*(njp1-1)
              ijk = ijkn - jo
              ijk1 = ijk - jo
              ijk2 = ik + jo*(2-1)
c              dfdyn(iblk,i,k) = (sc1(ijk) - sc1(ijk1)) / djv(i,nj)
              dfdyn(iblk,i,k) = (sc1(ijk) - sc1(ijk1)) / djv(ij1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 gamsc1(ijkn) = gamsc1(ijk2) 
              else if (bl_bc(iside,iblk) .eq. BCT_WALL) then
                 gamsc1(ijkn) = -gamsc1(ijk)
                 gamlam(ijkn) = -gamlam(ijk)
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE) then
                 gamsc1(ijkn) = gamsc1(ijk)
              endif
           enddo
        enddo

        do j = 2, nj
           do i = 2, ni
              ij = i + jo*(j-1) + ps_ijk
              ij1 = i + jo*(j-1) + pm_ij
              iside = 5
              ijkl = ij + ko*(1-1)
              ijk = ij + ko*(2-1)
              ijk1 = ijk + ko
              ijk2 = ij + ko*(nk-1)
              dfdzl(iblk,i,j) = (sc1(ijk1) - sc1(ijk)) /
     &                          (rp(ij1)*dkp(2))
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 gamsc1(ijkl) = gamsc1(ijk2) 
              else if (bl_bc(iside,iblk) .eq. BCT_WALL) then
                 gamsc1(ijkl) = -gamsc1(ijk)
                 gamlam(ijkl) = -gamlam(ijk)
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE) then
                 gamsc1(ijkl) = gamsc1(ijk)
              endif
              iside = 6
              ijkr = ij + ko*(nkp1-1)
              ijk = ijkr - ko
              ijk1 = ijk - ko
              ijk2 = ij + ko*(2-1)
              dfdzr(iblk,i,j) = (sc1(ijk) - sc1(ijk1)) /
     &                          (rp(ij1)*dkp(nk-1))
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 gamsc1(ijkr) = gamsc1(ijk2) 
              else if (bl_bc(iside,iblk) .eq. BCT_WALL) then
                 gamsc1(ijkr) = -gamsc1(ijk)
                 gamlam(ijkr) = -gamlam(ijk)
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE) then
                 gamsc1(ijkr) = gamsc1(ijk)
              endif
           enddo
        enddo

 10   continue

      return
      end
