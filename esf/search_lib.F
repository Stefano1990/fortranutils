C=======================================================================
C Update density field
C=======================================================================

      subroutine search_lib(sc1, sc1var, rhonp1, Temperature,
     &                    work,offwork,
     &                    istep, alpha, rho,rholo,rhohi)

      implicit none

C PARAMETERS
#include "dim.inc"
#include "cb.lgrid.inc"

      real     sc1(*), rhonp1(*), sc1var(*), Temperature(*)
      integer  offwork(*)
      real     work(*)
      integer  istep
      real     alpha,rholo,rhohi
      real     rho(*)
C LOCAL

      integer iblk
      integer nip1, njp1, nkp1 
      integer p3   

C-----------------------------------------------------------------------
C Update density
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)


C POINTERS
        p3  = mg_off(iblk,1) + 1
        call search_lib_block(sc1(p3),sc1var(p3),rhonp1(p3),
     &                    Temperature(p3),
     &                    nip1, njp1, nkp1, istep,
     &                    alpha, rho(p3),rholo,rhohi)
 10   continue

C-----------------------------------------------------------------------
C PERIODIC BOUNDARIES AND INTERFACES TREAMENTS
C-----------------------------------------------------------------------

C 'NOT' Necessary, because they are all definite values of mixture fraction 

C  Sc1 in each cell. Since periodic and/or interface treatment has already

C  been carried out for Sc1, it is not useful to do such treatment for 

C  rhos, T, and Ys. 

C  This also means 'NO' mpi operations in parallel version.

      call swap_sc(0,1,sc1,work,2)

      call swap_sc(0,1,rhonp1,work,2)
      call swap_sc(0,1,Temperature,work,2)

      return
      end
