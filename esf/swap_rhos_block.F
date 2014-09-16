C=======================================================================
C 1. rho -> rhonm1
C 2. rhonp1 -> rho
C=======================================================================

      subroutine swap_rhos_block(rho,rhonp1,rhonm1,nip1, njp1, nkp1)

      implicit none

C Input/Output Parameters

      integer  nip1, njp1, nkp1
      real     rho(nip1*njp1*nkp1)
      real     rhonp1(nip1*njp1*nkp1)
      real     rhonm1(nip1*njp1*nkp1)

C LOCAL

      integer  jo, ko, ni, nj, nk
      integer  i, j, k, jk, ijk
      integer  koff

C SOME CONSTANTS

      jo = nip1
      ko = nip1*njp1
      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C-----------------------------------------------------------------------
C LOOP OVER INTERIOR POINTS OF A BLOCK
C-----------------------------------------------------------------------

      do k = 1,nkp1
        koff = ko*(k-1)
        do j = 1,njp1
          jk = jo*(j-1) + koff
          do i = 1,nip1

            ijk   = i + jk

            rhonm1(ijk) = rho(ijk)
            rho(ijk) = rhonp1(ijk)

          end do
        end do
      end do

      return
      end
