C=======================================================================
C  SUBTRACT THE PRESSURE VALUE AT THE REFERENCE POINT FROM THE PRESSURE
C  FIELD
C=======================================================================

      subroutine p_ref_subtract(p)

      implicit none

#include "mpif.h"

#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

C INPUT/OUTPUT PARAMETERS

      real    p(*)

C LOCAL QUANTITIES

      integer iblk
      integer nip1, njp1, nkp1
      integer p3
      real    p_ref
      integer ierr

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C HERE HARD-CODED:
!  This needs to be generalised for a consistent point in the domain
!  At present it takes the 2,2,2 point of the first block on the first
! process.
      bl_ref_block = 65 
      bl_ref_ijk(1) = 2
      bl_ref_ijk(2) = 2
      bl_ref_ijk(3) = 2

C GET THE VALUE AT THE REFERENCE POINT

      if(pe.eq.bl2pe(bl_ref_block))then
        nip1 = mg_nijkp1(1,bl_ref_block,1)
        njp1 = mg_nijkp1(2,bl_ref_block,1)
        p3    = mg_off(bl_ref_block,1)
        p_ref = p(p3+bl_ref_ijk(1)+nip1*(bl_ref_ijk(2)-1)+
     &          nip1*njp1*(bl_ref_ijk(3)-1))
       if(debugging) then
       write(6,*) pe, iblkl2g(bl_ref_block),'reference pressure: ',p_ref
       endif
      endif
 
      CALL MPI_BCAST(p_ref, 1, REAL_TYPE_MPI, bl2pe(bl_ref_block),
     &MPI_COMM_WORLD, ierr)


C SUBTRACT THE REFERENCE VALUE

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

        p3 = mg_off(iblk,1) + 1

        call subtract_p(p(p3),p_ref,nip1,njp1,nkp1)

 10   continue

      return
      end

C-----------------------------------------------------------------------

      subroutine subtract_p(p,p_ref,nip1,njp1,nkp1)

      implicit none

      integer i, j, k
      integer nip1, njp1, nkp1
      real    p(nip1,njp1,nkp1)
      real    p_ref

      do k = 1,nkp1
        do j = 1,njp1
          do i = 1,nip1

            p(i,j,k) = p(i,j,k) - p_ref

          end do
        end do
      end do

      return
      end
