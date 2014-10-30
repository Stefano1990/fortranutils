C=======================================================================
C APPLY THE CONVECTIVE BOUNDARY CONDITION TO THE EXIT BOUNDARY TO OBTAIN
C THE VELOCITY AT THE NEW TIME STEP
C=======================================================================

      subroutine conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &  istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &  u,u_int,rho,dip,rp,dkp,
     &  dtime,uave)

      implicit none

#include "les.inc"

      integer d1, iside
      integer nip1, njp1, nkp1
      integer istr, iend, jstr, jend, kstr, kend
      integer io, jo, ko
      real    u(nip1*njp1*nkp1), u_int(nip1*njp1*nkp1),
     &        rho(nip1*njp1*nkp1)
      real    dip(nip1*njp1), rp(nip1*njp1), dkp(nkp1)
      real    dtime, uave

C LOCAL

      integer i, j, k, ijk, ijk2, ij, i1, k1
      integer ioff2, ioff3, joff, koff, off_ijk, off_ijk2
      real    const1, dx_inv
      real    rho1, rho2
      
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      const1 = uave*dtime

      if (d1 .ne. SV_W) then

        if (mod(iside,2) .eq. 1) then
          i  = istr - 1
          i1 = istr
          off_ijk  = -io
          off_ijk2 = 0
          const1   = -const1
        else
          i  = iend + 1
          i1 = iend
          off_ijk  = io
          off_ijk2 = io
        end if
        ioff3 = 1 + io*(i -1)
        ioff2 = 1 + io*(i1-1)

        do k = kstr,kend
          koff = ioff3 + ko*(k-1)
          do j = jstr,jend
            ijk  = koff  + jo*(j-1)
            ijk2 = ijk - off_ijk2
            ij   = ioff2 + jo*(j-1)
            rho1 = 0.5*(rho(ijk2)+rho(ijk2+off_ijk))
            rho2 = 0.5*(rho(ijk2)+rho(ijk2-off_ijk))
            dx_inv = 1.0/dip(ij)
cwk            u(ijk) =
cwk     &        (u_int(ijk)*rho1 + const1*dx_inv*u(ijk-off_ijk)*rho2)/
cwk     &        ((1.0 + const1*dx_inv)*rho1)
            u(ijk) = u_int(ijk)
          end do
        end do

      else

        if (mod(iside,2) .eq. 1) then
          k  = kstr - 1
          k1 = kstr
          off_ijk  = -ko
          off_ijk2 = 0
          const1   = -const1
        else
          k  = kend + 1
          k1 = kend
          off_ijk  = ko
          off_ijk2 = ko
        end if
        koff = ko*(k-1)

        do j = jstr,jend
          joff = koff + jo*(j-1)
          do i = istr,iend
            ijk  = i + joff
            ijk2 = ijk - off_ijk2
            ij   = i + nip1*(j-1)
            dx_inv = 1.0/(rp(ij)*dkp(k1))
            rho1 = 0.5*(rho(ijk2)+rho(ijk2+off_ijk))
            rho2 = 0.5*(rho(ijk2)+rho(ijk2-off_ijk))
cwk            u(ijk) =
cwk     &        (u_int(ijk)*rho1 + const1*dx_inv*u(ijk-off_ijk)*rho2)/
cwk     &        ((1.0 + const1*dx_inv)*rho1)
            u(ijk) = u_int(ijk)
          end do
        end do

      end if

      return
      end
