C=======================================================================
C CALCULATE MASS FLUX FOR A GIVEN PLANE
C=======================================================================

C NOTE:
C 1. FOR CURVILINEAR ORTHORGONAL GRID USED BY ASALES
C 2. MUST MAINTAIN CONSISTENCY WITH THE WAY OF CALCULATING mass_flux
C    USED IN SUBROUTINE sour_p_block
C 3. VIRTUAL DENSITY MAY BE USED IF THE PLANE IS A BOUNDARY. BUT NO
C    SUCH DENSITY HAS BEEN DEFINED AT THIS MOMENT.

      subroutine mass_flux(plane,nip1,njp1,nkp1,d1,iva,
     &  istr,iend,jstr,jend,kstr,kend,
     &  area,area_rho,mflux,u,v,w,rho,
     &  rdju,rdiv,didjp,dkp)

      implicit none

#include "les.inc"

      integer plane, d1
      integer nip1, njp1, nkp1
      integer iva
      integer istr, iend, jstr, jend, kstr, kend
      real    u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real    rho(nip1*njp1*nkp1)
      real    rdju(nip1*njp1), rdiv(nip1*njp1), didjp(nip1*njp1)
      real    area, mflux,area_rho
      real    dkp(nkp1)

C LOCAL

      integer jo, ko, joff, koff, ikoff, jkoff
      integer i, j, k, ijk
      real    areaf

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C INTIAILIZE

      area  = 0.
      mflux = 0.
      area_rho = 0.0

C NOT THE VARIABLE CONSIDERED

      if (d1 .ne. iva) return

C-----------------------------------------------------------------------
C MASS FLUX AND AREA
C-----------------------------------------------------------------------

      jo = nip1
      ko = nip1*njp1

      if (d1 .eq. SV_U) then

        do k = kstr,kend
          ikoff = ko*(k-1) + plane
          do j = jstr,jend
            ijk   = jo*(j-1) + ikoff
            areaf = rdju(plane+jo*(j-1))*dkp(k)
            area  = area  + areaf
            mflux = mflux + 0.5*(rho(ijk)+rho(ijk-1))*u(ijk)*areaf
            area_rho = area_rho + 0.5*(rho(ijk)+rho(ijk-1))*areaf
          end do
        end do

      else if (d1 .eq. SV_V) then

        joff = jo*(plane-1)
        do k = kstr,kend
          jkoff = ko*(k-1) + joff
          do i = istr,iend
            ijk   = i + jkoff
            areaf = rdiv(i+joff)*dkp(k)
            area  = area  + areaf
            mflux = mflux + 0.5*(rho(ijk)+rho(ijk-jo))*v(ijk)*areaf
            area_rho = area_rho + 0.5*(rho(ijk)+rho(ijk-jo))*areaf
          end do
        end do

      else

        koff = ko*(plane-1)
        do j = jstr,jend
          joff  = jo*(j-1)
          jkoff = joff + koff
          do i = istr,iend
            ijk   = i + jkoff
            areaf = didjp(i+joff)
            area  = area + areaf
            mflux = mflux + 0.5*(rho(ijk)+rho(ijk-ko))*w(ijk)*areaf
            area_rho = area_rho + 0.5*(rho(ijk)+rho(ijk-ko))*areaf
          end do
        end do

      end if

      return
      end
