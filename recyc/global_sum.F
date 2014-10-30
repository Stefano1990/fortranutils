C=======================================================================
C CALCULATE GLOBAL SUM OF A GRID QUANTITY FOR ALL BLOCKS
C=======================================================================

      subroutine global_sum(ilev,mg_or_sg,f,sum,termmax,termmin)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "const.inc"
#include "comm.inc"
#include "data_type.inc"

C INCLUDE
#include "cb.lgrid.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"

C INPUT/OUTPUT PARAMETERS

      integer ilev, mg_or_sg
      real    f(*)
      real    sum, termmax, termmin

C LOCAL QUANTITIES

      integer iblk, iblk0, ierr
      integer nip1, njp1, nkp1
      integer is, ie, js, je, ks, ke
      real    sum_blk, sum2
      integer p3, p3x
      integer status(MPI_STATUS_SIZE), i
#ifdef DEBUG
      real sum3(maxbl)
#endif

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

#ifdef DEBUG
      do iblk=1,numbl
        sum3(iblk) = 0.0
      end do
#endif
 
      sum = 0.
      termmax = -large_number
      termmin =  large_number

      do 10 iblk = 1,numbl

        p3 = mg_off(iblk,ilev) + 1
        if (mg_or_sg .eq. 1) then
          p3x = p3
        else
          p3x = sg_off(iblk) + 1
        end if

        nip1 = mg_nijkp1(1,iblk,ilev)
        njp1 = mg_nijkp1(2,iblk,ilev)
        nkp1 = mg_nijkp1(3,iblk,ilev)

        is = nijkse(1,iblk,ilev)
        ie = nijkse(2,iblk,ilev)
        js = nijkse(3,iblk,ilev)
        je = nijkse(4,iblk,ilev)
        ks = nijkse(5,iblk,ilev)
        ke = nijkse(6,iblk,ilev)

        call global_sum_block(f(p3),nip1,njp1,nkp1,
     &       is,ie,js,je,ks,ke,sum_blk,termmax,termmin)

        sum = sum + sum_blk

#ifdef DEBUG
        sum3(iblk) = sum3(iblk) + sum_blk
#endif

 10   continue

#ifdef DEBUG
      call usr_sum_double (sum3, sum, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
      sum2 = sum
      call MPI_ALLREDUCE(sum2, sum, 1, REAL_TYPE_MPI,
     &MPI_SUM, MPI_COMM_WORLD, ierr)
#endif



      return
      end



C=======================================================================
C CALCULATE NORMS FOR A BLOCK
C=======================================================================

      subroutine global_sum_block(f,nip1,njp1,nkp1,
     &  nis,nie,njs,nje,nks,nke,sum,termmax,termmin)

      implicit none

C INPUT/OUTPUT PARAMETERS

      integer nip1, njp1, nkp1
      integer nis, nie, njs, nje, nks, nke
      real    f(nip1*njp1*nkp1)
      real    sum, termmax, termmin

C LOCAL QUANTITIES

      integer i, j, k, ijk
      integer io, jo, ko, koff, jkoff

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      sum = 0.

      io   = 1
      jo   = nip1
      ko   = nip1*njp1

      do k = nks,nke
        koff = ko*(k-1)
        do j = njs,nje
          jkoff = jo*(j-1) + koff
          do i = nis,nie

            ijk = i + jkoff
            sum = sum + f(ijk)
            termmax = max(termmax,f(ijk))
            termmin = min(termmin,f(ijk))

          end do
        end do
      end do

      return
      end
