C=======================================================================
C CALCULATE scalar(mixture fraction) variance
C=======================================================================

      subroutine set_sc1var_block(u,v,w,sc1,rho,
     &             dip,djp,dkp,dkw,rp,sc1var,sc1dis)

      implicit none

C PARAMETERS

#        include "dim.inc"

C INCLUDE FILES

#        include "cb.gridp.inc"
#        include "cb.lgridp.inc"
#        include "cb.lgrid.inc"
#        include "cb.control.inc"

C Input/Output Parameters

      real     rho(nip1*njp1*nkp1)
      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     sc1(nip1*njp1*nkp1)
      real     sc1var(nip1*njp1*nkp1),sc1dis(nip1*njp1*nkp1)
      real     dip(nip1,njp1),djp(nip1,njp1)
      real     dkp(nkp1), dkw(nkp1)
      real     rp(nip1,njp1)

C LOCAL

      integer  i, j, k, jk, ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr
      integer  koff
      real     di, dj, dk
      real     radp
      real     delta_xyz
      real     dfdip, dfdjp, dfdkp, dfdijk
      real     fwf, fef, fsf, fnf, flf, frf
      real     cfv

C-----------------------------------------------------------------------
      data     cfv/0.1/
C-----------------------------------------------------------------------

C SOME CONSTANTS

      jo = nip1
      ko = nip1*njp1
      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C-----------------------------------------------------------------------
C LOOP OVER INTERIOR POINTS OF A BLOCK
C-----------------------------------------------------------------------

      do k = 2,nk
        koff = ko*(k-1)
        do j = 2,nj
          jk = jo*(j-1) + koff
          do i = 2,ni

            ijk   = i + jk

            ijkw  = ijk  - 1
            ijke  = ijk  + 1
            ijks  = ijk  - jo
            ijkn  = ijk  + jo
            ijkl  = ijk  - ko
            ijkr  = ijk  + ko

            fwf = 0.5*(sc1(ijk)+sc1(ijkw))
            fef = 0.5*(sc1(ijk)+sc1(ijke))
            fsf = 0.5*(sc1(ijk)+sc1(ijks))
            fnf = 0.5*(sc1(ijk)+sc1(ijkn))
            flf = 0.5*(sc1(ijk)+sc1(ijkl))
            frf = 0.5*(sc1(ijk)+sc1(ijkr))

            dfdip = (fef-fwf)/dip(i,j)
            dfdjp = (fnf-fsf)/djp(i,j)
            dfdkp = (frf-flf)/(rp(i,j)*dkw(k))

            dfdijk = dfdip**2 + dfdjp**2 + dfdkp**2

C+++++++++++++++++++++++++++++Important++++++++++++++++++++++++++++
            di = dip(i,j)
            dj = djp(i,j)
            radp = rp(i  ,j  )
            dk = radp*dkp(k)

            delta_xyz = (di*dj*dk)**.33333333
C+++++++++++++++++++++++++++++Important++++++++++++++++++++++++++++

            sc1var(ijk) = cfv*dfdijk*(delta_xyz**2)


          end do
        end do
      end do


      return
      end
