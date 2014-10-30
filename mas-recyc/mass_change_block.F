      subroutine mass_change_block(rho,rhonp1,mass_chg_blk,
     &                             rdidjp,dkw, dtime)

      implicit none

C PARAMETERS

#include "dim.inc"

C INCLUDE FILES

#include "cb.gridp.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"

C Input/Output Parameters

      real     rho(nip1*njp1*nkp1)
      real     rhonp1(nip1*njp1*nkp1)
      real     mass_chg_blk

      real     rdidjp(nip1,njp1) 
      real     dkw(nkp1)
      real     dtime, ddtime
C LOCAL

      integer  i, j, k, jk, ijk
      integer  koff
      real     vol 
      real     cell_change

C SOME CONSTANTS

      jo = nip1
      ko = nip1*njp1
      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------
       ddtime = 1.0 / dtime
       mass_chg_blk = 0.0

C-----------------------------------------------------------------------
C LOOP OVER INTERIOR POINTS OF A BLOCK
C-----------------------------------------------------------------------

CX      do k = 1,nkp1
      do k = 2,nk
        koff = ko*(k-1)
CX        do j = 1,njp1
        do j = 2,nj
          jk = jo*(j-1) + koff
CX          do i = 1,nip1
          do i = 2,ni
            ijk   = i + jk

            vol = rdidjp(i,j)*dkw(k)
            cell_change = ddtime * vol*(rhonp1(ijk)-rho(ijk))
            mass_chg_blk = mass_chg_blk + cell_change

          end do
        end do
      end do

      return
      end
