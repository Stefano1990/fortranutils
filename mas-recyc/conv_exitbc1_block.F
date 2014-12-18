C=======================================================================
C APPLY THE CONVECTIVE BOUNDARY CONDITION TO THE EXIT BOUNDARY TO OBTAIN
C THE INTERMEDIATE VELOCITY
C=======================================================================

      subroutine conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &  istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &  u,u_int,rho,rhonp1,hn,dip,rp,dkp,
     &  dtime,uave)

      implicit none

#include "les.inc"

      integer d1, iside
      integer nip1, njp1, nkp1
      integer istr, iend, jstr, jend, kstr, kend
      integer io, jo, ko
      real    u(nip1*njp1*nkp1), u_int(nip1*njp1*nkp1),
     &        hn(nip1*njp1*nkp1), rho(nip1*njp1*nkp1)
      real    dip(nip1*njp1), rp(nip1*njp1), dkp(nkp1)
      real    dtime, uave

C LOCAL

      integer i, j, k, ijk, ijk2, ij, i1, k1
      integer ioff2, ioff3, joff, koff, off_ijk, off_ijk2
      real    const1, const2, h, uave1
      real    rho1, rho2


      real    rhonp1(nip1*njp1*nkp1)
      real    rhonp1_1


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      const1 = 1.5*dtime
      const2 = 0.5*dtime

C-----------------------------------------------------------------------
C FOR THE u-, v- VELOCITIES
C-----------------------------------------------------------------------

      if (d1 .ne. SV_W) then

        if (mod(iside,2) .eq. 1) then
          i  = istr - 1
          i1 = istr
          off_ijk  = -io
          off_ijk2 = 0
          uave1 = -uave
        else
          i  = iend + 1
          i1 = iend
          off_ijk  = io
          off_ijk2 = io
          uave1 = uave
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
            h   = uave1*(rho1*u(ijk)-rho2*u(ijk-off_ijk)) / dip(ij)

            rhonp1_1 = 0.5*(rhonp1(ijk2)+rhonp1(ijk2+off_ijk))
       u_int(ijk) = (rho1*u(ijk) - const1*h + const2*hn(ijk))/rhonp1_1

!these comments are for a neumann condition
!should really formalise this with a switch

!       u_int(ijk) = u(ijk-1)  
!this can help with negative velocities at the outflow boundary
!       if(d1 .eq. SV_U .and. u_int(ijk) .lt. 0.0) u_int(ijk) = 0.0

            hn(ijk) = h
          end do
        end do

C-----------------------------------------------------------------------
C FOR THE w-VELOCITY
C-----------------------------------------------------------------------

      else

        if (mod(iside,2) .eq. 1) then
          k  = kstr - 1
          k1 = kstr
          off_ijk  = -ko
          off_ijk2 = 0
          uave1 = -uave
        else
          k  = kend + 1
          k1 = kend
          off_ijk  = ko
          off_ijk2 = ko
          uave1 = uave
        end if
        koff = ko*(k-1)

        do j = jstr,jend
          joff = koff + jo*(j-1)
          do i = istr,iend
            ijk  = i + joff
            ijk2 = ijk - off_ijk2
            ij   = i + nip1*(j-1)
            rho1 = 0.5*(rho(ijk2)+rho(ijk2+off_ijk))
            rho2 = 0.5*(rho(ijk2)+rho(ijk2-off_ijk))
            h    = uave1*(rho1*u(ijk)-rho2*u(ijk-off_ijk))
     &           /(rp(ij)*dkp(k1))

            rhonp1_1 = 0.5*(rhonp1(ijk2)+rhonp1(ijk2+off_ijk))
       u_int(ijk) = (rho1*u(ijk) - const1*h + const2*hn(ijk))/rhonp1_1


!these comments are for a neumann condition
!should really formalise this with a switch       
!       u_int(ijk) = u(ijk - 1)

            hn(ijk) = h
          end do
        end do

      end if

      return
      end
