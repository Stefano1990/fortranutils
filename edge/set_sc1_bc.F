C=======================================================================
C LOOP THROUGH ALL BLOCKS TO SET UP BOUNDARY VALUES FOR THE SCALAR
C=======================================================================

      subroutine set_sc1_bc(sc1)

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

      real         sc1(*)
      real         sc1tmp
      
C LOCAL

      integer      iblk, iside, i, j, k, ij, ik, jk, jk1
      integer      nip1, njp1, nkp1, ni, nj, nk
      integer      io, jo, ko
      integer      ps_ijk, ps_ijk_last
      integer      ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr, ijk2
      integer      ijk3

C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        ps_ijk = sg_off(iblk)

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
              ijkw = 1 + jk
              ijk = 2 + jk
              ijk2 = ni + jk
              ijk3 = nip1 + jk
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 sc1(ijkw) = sc1(ijk2) 
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE .and.
     &            bl_bc(iside,iblk) .ne. BCT_FIXED) then
                 sc1(ijkw) = sc1(ijk)
              endif
C.. TEMPORARY
c              if (iblk .eq. 1) then
c                 sc1(ijkw) = 2.0 - sc1(ijk)
c              else if (iblk .eq. 2) then
c                 sc1(ijkw) = -sc1(ijk)
c              endif

              iside = 2
              ijke = nip1 + jk
              ijk = ijke - 1
              ijk2 = 2 + jk
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 sc1(ijke) = sc1(ijk2) 
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE .and.
     &            bl_bc(iside,iblk) .ne. BCT_FIXED) then
                 sc1(ijke) = sc1(ijk)
              endif
           enddo
        enddo

        do i = 2, ni
           do k = 2, nk
              ik = i + ko*(k-1) + ps_ijk 
              iside = 3
              ijks = ik + jo*(1-1)
              ijk = ik + jo*(2-1)
              ijk2 = ik + jo*(nj-1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 sc1(ijks) = sc1(ijk2) 
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE .and.
     &            bl_bc(iside,iblk) .ne. BCT_FIXED) then
                 sc1(ijks) = sc1(ijk)
              endif
C.. TEMPORARY
c              if (iblk .eq. 2) then
c                 sc1(ijks) = -sc1(ijk)
c              endif
              iside = 4
              ijkn = ik + jo*(njp1-1)
              ijk = ijkn - jo
              ijk2 = ik + jo*(2-1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 sc1(ijkn) = sc1(ijk2) 
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE .and.
     &            bl_bc(iside,iblk) .ne. BCT_FIXED) then
                 sc1(ijkn) = sc1(ijk)
              endif
C.. TEMPORARY
c              if (iblk .eq. 1) then
c                 sc1(ijkn) = 2.0 - sc1(ijk)
c              else if (iblk .eq. 2) then
c                 sc1(ijkn) = -sc1(ijk)
c              endif
           enddo
        enddo

        do j = 2, nj
           do i = 2, ni
              ij = i + jo*(j-1) + ps_ijk
              iside = 5
              ijkl = ij + ko*(1-1)
              ijk = ij + ko*(2-1)
              ijk2 = ij + ko*(nk-1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 sc1(ijkl) = sc1(ijk2) 
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE .and.
     &            bl_bc(iside,iblk) .ne. BCT_FIXED) then
                 sc1(ijkl) = sc1(ijk)
              endif
              iside = 6
              ijkr = ij + ko*(nkp1-1)
              ijk = ijkr - ko
              ijk2 = ij + ko*(2-1)
              if (bl_bc(iside,iblk) .eq. BCT_PERIODIC) then
                 sc1(ijkr) = sc1(ijk2) 
              else if (bl_bc(iside,iblk) .ne. BCT_INTERFACE .and.
     &            bl_bc(iside,iblk) .ne. BCT_FIXED) then
                 sc1(ijkr) = sc1(ijk)
              endif
           enddo
        enddo

 10   continue

      return
      end
