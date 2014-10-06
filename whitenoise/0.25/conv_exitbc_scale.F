C=======================================================================
C SCALE THE VELOCITY AT EXIT BOUNDARIES
C=======================================================================

      subroutine conv_exitbc_scale(iva,iside,iblk,nip1,njp1,nkp1,
     &  u,scalef)

      implicit none

#include "dim.inc"

#include "cb.lcontrol.inc"

      integer iva, iside, iblk
      integer nip1, njp1, nkp1
      real    u(nip1*njp1*nkp1)
      real    scalef
      
C LOCAL

      integer d1, d2, d3
      integer io, jo, ko, ioff, koff
      integer i, j, k, ijk
      integer istr, iend, jstr, jend, kstr, kend
      integer ijko(3), ijkse(6)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      d1 = (iside+1)/2
      d2 = cyclic1(1,d1)
      d3 = cyclic1(2,d1)

      ijko(1) = 1
      ijko(2) = nip1
      ijko(3) = nip1*njp1

      io = ijko(d1)
      jo = ijko(d2)
      ko = ijko(d3)

      call set_loop_bound(iva,1,1,iblk,iblk,1,ijkse)

      istr = ijkse(2*d1-1)
      iend = ijkse(2*d1  )
      jstr = ijkse(2*d2-1)
      jend = ijkse(2*d2  )
      kstr = ijkse(2*d3-1)
      kend = ijkse(2*d3  )

      if (mod(iside,2) .eq. 1) then
        i = istr - 1
      else
        i = iend + 1
      end if

      ioff = 1 + io*(i-1)

      do k = kstr,kend
        koff = ioff + ko*(k-1)
        do j = jstr,jend
          ijk = koff + jo*(j-1)
          u(ijk) = u(ijk)*scalef
C  Temporary
c          u(ijk) = amax1(0.0, u(ijk))
        end do
      end do

      return
      end
