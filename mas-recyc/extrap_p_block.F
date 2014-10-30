C=======================================================================
C  EXTRAPOLATE THE PRESSURE TO THE BOUNDARY POINTS
C=======================================================================

      subroutine extrap_p_block(p,nip1,njp1,nkp1,bc)


      implicit none

#include "bc_type.inc"

C INPUT/OUTPUT PARAMETERS

      integer nip1, njp1, nkp1
      integer bc(6)
      real    p(nip1*njp1*nkp1)

C LOCAL QUANTITIES

      integer d1, d2, d3
      integer iside, i, j, k, jk, ijk, io, io1, jo, ko, jko
      integer ijko(3), ijkbc(6), nijkp1(3), co(2,3)
      integer ioff, koff, jkoff

      data    co/2,3,1,3,1,2/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      ijko(1)   = 1
      ijko(2)   = nip1
      ijko(3)   = nip1*njp1
      ijkbc(1)  = 1
      ijkbc(2)  = nip1
      ijkbc(3)  = 1
      ijkbc(4)  = njp1
      ijkbc(5)  = 1
      ijkbc(6)  = nkp1
      nijkp1(1) = nip1
      nijkp1(2) = njp1
      nijkp1(3) = nkp1

      do 10 iside = 1,6

        if (bc(iside).ne.BCT_INTERFACE .and.
     &      bc(iside).ne.BCT_PERIODIC) then

          d1 = (iside+1)/2
          d2 = co(1,d1)
          d3 = co(2,d1)

          io = ijko(d1)
          if(mod(iside,2).eq.0) then
            io1 = -io
          else
            io1 =  io
          end if
          jo = ijko(d2)
          ko = ijko(d3)

          i = ijkbc(iside)

          ioff = 1 + io*(i-1)
          do k = 1,nijkp1(d3)
            koff = ioff + ko*(k-1)
            do j = 1,nijkp1(d2)
              ijk = jo*(j-1) + koff
              p(ijk) = p(ijk+io1)
            end do
          end do

        end if

 10   continue

      return
      end
