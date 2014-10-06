C=======================================================================
C INITIALIZE AN ARRAY TO A GIVEN VALUE
C=======================================================================

      subroutine init_array(f,value,mg_or_sg,lev1,lev2)

      implicit none

C PARAMETERS

#include "dim.inc"

C INPUT/OUTPUT

      integer mg_or_sg, lev1, lev2
      real    f(*), value

C LOCAL

      integer i, is, ie, ilev, iblk
      integer nip1, njp1, nkp1
      integer npts, poff

C INCLUDE

#include "cb.lgrid.inc"
#include "cb.control.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      if (mg_or_sg.eq.1) then

        do ilev = lev1,lev2
          do iblk = 1,numbl

            nip1 = mg_nijkp1(1,iblk,ilev)
            njp1 = mg_nijkp1(2,iblk,ilev)
            nkp1 = mg_nijkp1(3,iblk,ilev)

            npts = nip1*njp1*nkp1

            poff = mg_off(iblk,ilev)
            is = poff + 1
            ie = poff + npts

            do i = is,ie
              f(i) = value
            end do

          end do
        end do

      else

        do iblk = 1,numbl

          nip1 = bl_nijkp1(1,iblk)
          njp1 = bl_nijkp1(2,iblk)
          nkp1 = bl_nijkp1(3,iblk)

          npts = nip1*njp1*nkp1

          poff = sg_off(iblk)
          is = poff + 1
          ie = poff + npts

          do i = is,ie
            f(i) = value
          end do

        end do

      end if

      return
      end
