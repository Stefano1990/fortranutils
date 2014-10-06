C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+ Update density and scalar fields
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine search_lib_block(sc1, sc1var, rhonp1, Temperature,
     &                    nip1, njp1, nkp1, istep, alpha,rho,
     &                    rholo,rhohi)

      implicit none

C ARGUMENTS

      integer nip1, njp1, nkp1, istep
      real  sc1(nip1*njp1*nkp1), sc1var(nip1*njp1*nkp1)
      real  rhonp1(nip1*njp1*nkp1), Temperature(nip1*njp1*nkp1)
	  real   alpha,rholo,rhohi
      real  rho(nip1*njp1*nkp1)

C Local
      integer  i, j, k, jk, ijk
      integer  koff
      integer  Zi, Zj
      integer  indexj, indexi
      integer  jo, ko, ni, nj, nk   
      real     rho_tmp, T_tmp, rho_mix

C SOME CONSTANTS

      jo = nip1
      ko = nip1*njp1
      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1


      do k = 1,nkp1
        koff = ko*(k-1)
        do j = 1,njp1
          jk = jo*(j-1) + koff
          do i = 1,nip1

            ijk   = i + jk

              rho_mix=1.0/(sc1(ijk)/rhohi + (1.0-sc1(ijk))/rholo)
              if (sc1(ijk) .gt. 1.0) rho_mix = rhohi 
              if (sc1(ijk) .lt. 0.0) rho_mix = rholo  

              rho_tmp = rho_mix

C perform relaxation operation
      rhonp1(ijk) = alpha*rho(ijk) + (1.0 -alpha)*rho_tmp
C update scalar field
      sc1(ijk) = (sc1(ijk)*rho(ijk))/rhonp1(ijk)

C set limits on scalar field
      if (sc1(ijk) .gt. 1.0) sc1(ijk) =1.0


          end do
        end do
      end do

      return
      end

