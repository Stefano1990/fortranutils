C=======================================================================
C CORRECT THE INTERMEDIATE VELOCITY (ADAMS-BASHFORTH SCHEME) FOR A BLOCK
C=======================================================================

      subroutine uvw_adv_block(ivar,nip1,njp1,nkp1,istr,jstr,kstr,
     &  dtime,diu,rp,dkw,u,u_int,rho,p,blk,step)

      implicit none

#include "les.inc"

C INPUT/OUTPUT PARAMETERS

      integer ivar,blk,step
      integer nip1, njp1, nkp1
      integer istr, jstr, kstr
      real    dtime
      real    diu(nip1*njp1), rp(nip1*njp1)
      real    dkw(nkp1)
      real    u(nip1*njp1*nkp1), u_int(nip1*njp1*nkp1),
     &        rho(nip1*njp1*nkp1), p(nip1*njp1*nkp1)

C LOCAL QUANTITIES
 
      integer i, j, k
      integer ijk, ijkw, ijkl, ij
      integer koff, jkoff, joff_ij
      integer ni, nj, nk
      integer io, jo, ko, io_or_jo
      real    ddtime
      real    dpdi, rhou

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

      io = 1
      jo = nip1
      ko = nip1*njp1

      ddtime = 1.5*dtime

C FOR THE u- OR v-VELOCITY COMPONENT

      if(ivar .ne. SV_W) then

        if (ivar .eq. SV_U) then
          io_or_jo = io
        else
          io_or_jo = jo
        end if
!      if(ivar .eq. SV_V) then
!        write(*,*) blk,jstr,jo,io_or_jo,njp1-1
!      endif

        do k = kstr,nk
          koff = ko*(k-1)
          do j = jstr,nj
            jkoff   = jo  *(j-1) + koff
            joff_ij = nip1*(j-1)
            do i = istr,ni

              ijk  = i   + jkoff
              ijkw = ijk - io_or_jo

              ij = i + joff_ij

              dpdi = (p(ijk) - p(ijkw))/diu(ij)
              rhou = 0.5*(rho(ijk) + rho(ijkw))
             u(ijk) = u_int(ijk) - ddtime*dpdi/rhou
!        if (ivar .eq. SV_U) then
!          if(u(ijk) .lt. 0.0) then
!             u(ijk) = 4.0
!          endif
!        endif
            end do
          end do
        end do

      else

!        write(*,*) ivar,blk

        do k = kstr,nk
          koff = ko*(k-1)
          do j = jstr,nj
            jkoff   = jo  *(j-1) + koff
            joff_ij = nip1*(j-1)
            do i = jstr,ni

              ijk  = i   + jkoff
              ijkl = ijk - ko

              ij = i + joff_ij

              dpdi = (p(ijk) - p(ijkl))/(rp(ij)*dkw(k))
              rhou = 0.5*(rho(ijk) + rho(ijkl))
              u(ijk) = u_int(ijk) - ddtime*dpdi/rhou
            end do
          end do
        end do

      end if

      return
      end
