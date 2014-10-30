C=======================================================================
C SET UP OFFSETS
C=======================================================================

C GRID LEVEL INDEX: THE FINEST GRID LEVEL IS LABELLED AS 1, THE NEXT
C   COARSER GRID LEVEL AS 2, AND SO ON

      subroutine set_offset()

      implicit none

#include "mpif.h"
C PARAMETERS

#include "dim.inc"
#include "bc_type.inc"
#include "data_type.inc"
#include "comm.inc"

C LOCAL

      integer ired, ierr

      integer iblk, lev
      integer nip1, njp1, nkp1
      integer nijkp1_smallest
      integer nummgl_found
      logical coarser

      data    nijkp1_smallest/4/

C INCLUDE

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C OFFSET FOR SINGLE GRID QUANTITIES
C-----------------------------------------------------------------------

      sg_off_1d_i(1)  = 0
      sg_off_1d_j(1)  = 0
      sg_off_1d_k(1)  = 0
      sg_off_2d_ij(1) = 0
      sg_off_2d_jk(1) = 0
      sg_off_2d_ki(1) = 0
      sg_off(1)       = 0

      do 10 iblk = 1,numbl-1

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        sg_off_1d_i(iblk+1)  = sg_off_1d_i(iblk)  + nip1
        sg_off_1d_j(iblk+1)  = sg_off_1d_j(iblk)  + njp1
        sg_off_1d_k(iblk+1)  = sg_off_1d_k(iblk)  + nkp1
        sg_off_2d_ij(iblk+1) = sg_off_2d_ij(iblk) + nip1*njp1
        sg_off_2d_jk(iblk+1) = sg_off_2d_jk(iblk) + njp1*nkp1
        sg_off_2d_ki(iblk+1) = sg_off_2d_ki(iblk) + nkp1*nip1
        sg_off(iblk+1)       = sg_off(iblk)  + nip1*njp1*nkp1

 10   continue

C-----------------------------------------------------------------------
C OFFSET FOR MULTIGRID GRID QUANTITIES
C-----------------------------------------------------------------------

      maxijk1 = 0
      maxijk2 = 0
      maxijk3 = 0

      do 15 iblk = 1,numbl
        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)
        mg_nijkp1(1,iblk,1) = nip1
        mg_nijkp1(2,iblk,1) = njp1
        mg_nijkp1(3,iblk,1) = nkp1
        maxijk1 = max(maxijk1,nip1,njp1,nkp1)
        maxijk2 = max(maxijk2,nip1*njp1,njp1*nkp1,nkp1*nip1)
        maxijk3 = max(maxijk3,nip1*njp1*nkp1)
        mg_dim1(iblk) = nip1
        mg_dim2(iblk) = max(nip1*njp1,njp1*nkp1,nkp1*nip1)
        mg_dim3(iblk) = nip1*njp1*nkp1
 15   continue

C FIND THE LARGEST MULTIGRID LEVELS USED FOR EACH BLOCK

      do 20 lev = 1,maxmgl-1
        coarser = .true.
        do 25 iblk = 1,numbl

          if (mod(mg_nijkp1(1,iblk,lev),2).eq.0 .and.
     &        mod(mg_nijkp1(2,iblk,lev),2).eq.0 .and.
     &        mod(mg_nijkp1(3,iblk,lev),2).eq.0 .and. coarser) then
            nip1 = (mg_nijkp1(1,iblk,lev))/2 + 1
            njp1 = (mg_nijkp1(2,iblk,lev))/2 + 1
            nkp1 = (mg_nijkp1(3,iblk,lev))/2 + 1

            if (nip1.ge.nijkp1_smallest .and.
     &          njp1.ge.nijkp1_smallest .and.
     &          nkp1.ge.nijkp1_smallest) then

              mg_nijkp1(1,iblk,lev+1) = nip1
              mg_nijkp1(2,iblk,lev+1) = njp1
              mg_nijkp1(3,iblk,lev+1) = nkp1
              mg_dim1(iblk) = mg_dim1(iblk) + nip1
              mg_dim2(iblk) = mg_dim2(iblk) +
     &                        max(nip1*njp1,njp1*nkp1,nkp1*nip1)
              mg_dim3(iblk) = mg_dim3(iblk) + nip1*njp1*nkp1

            else
              coarser = .false.
            end if
          else
            coarser = .false.
          end if

 25     continue

        if (coarser) then
          nummgl_found = lev + 1
        else
          goto 1000
        end if

 20   continue

 1000 continue

      ired = nummgl_found
      call MPI_ALLREDUCE(ired, nummgl_found, 1, INTEGER_TYPE_MPI,
     &  MPI_MIN, MPI_COMM_WORLD, ierr)

      if (nummgl .eq. 0) nummgl = nummgl_found
      if (master) then
        if (nummgl .eq. 0) then
          write(*,'(''*** Info: Num. grid levels found: '',i2,
     &              '' and used: '',i2)') nummgl_found,nummgl
        else if (nummgl .le. nummgl_found) then
          write(*,'(''*** Info: Num. grid levels found: '',i2,
     &              '' and used: '',i2)') nummgl_found,nummgl
        else
          write(*,'(''*** Error: Num. grid levels found: '',i2,
     &              '' but required: '',i2)') nummgl_found,nummgl
          stop
        end if
      end if

      maxijk1_mg = 0
      maxijk2_mg = 0
      maxijk3_mg = 0
      do iblk = 1,numbl
        maxijk1_mg = max(maxijk1_mg,mg_dim1(iblk))
        maxijk2_mg = max(maxijk2_mg,mg_dim2(iblk))
        maxijk3_mg = max(maxijk3_mg,mg_dim3(iblk))
      end do

C CALCULATE THE OFFSETS

      mg_off_1d_i(1,1)  = 0
      mg_off_1d_j(1,1)  = 0
      mg_off_1d_k(1,1)  = 0
      mg_off_2d_ij(1,1) = 0
      mg_off_2d_jk(1,1) = 0
      mg_off_2d_ki(1,1) = 0
      mg_off(1,1)       = 0

      mg_off_1d_z(1,1)  = 0
      mg_off_2d_xy(1,1) = 0

      do 30 iblk = 1,numbl-1

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

        mg_off_1d_i(iblk+1,1)  = mg_off_1d_i(iblk,1)  + nip1
        mg_off_1d_j(iblk+1,1)  = mg_off_1d_j(iblk,1)  + njp1
        mg_off_1d_k(iblk+1,1)  = mg_off_1d_k(iblk,1)  + nkp1
        mg_off_2d_ij(iblk+1,1) = mg_off_2d_ij(iblk,1) + nip1*njp1
        mg_off_2d_jk(iblk+1,1) = mg_off_2d_jk(iblk,1) + njp1*nkp1
        mg_off_2d_ki(iblk+1,1) = mg_off_2d_ki(iblk,1) + nkp1*nip1
        mg_off(iblk+1,1)       = mg_off(iblk,1) +  nip1*njp1*nkp1
        
C ONLY FOR THE ARRAYS STORING THE COORDINATES

        mg_off_1d_z(iblk+1,1)  = mg_off_1d_z(iblk,1) + nkp1+1
        mg_off_2d_xy(iblk+1,1) = mg_off_2d_xy(iblk,1)+(nip1+1)*(njp1+1)

 30   continue

      do 35 lev = 2,nummgl

        nip1 = mg_nijkp1(1,numbl,lev-1)
        njp1 = mg_nijkp1(2,numbl,lev-1)
        nkp1 = mg_nijkp1(3,numbl,lev-1)

        mg_off_1d_i(1,lev)  = mg_off_1d_i(numbl,lev-1)  + nip1
        mg_off_1d_j(1,lev)  = mg_off_1d_j(numbl,lev-1)  + njp1
        mg_off_1d_k(1,lev)  = mg_off_1d_k(numbl,lev-1)  + nkp1
        mg_off_2d_ij(1,lev) = mg_off_2d_ij(numbl,lev-1) + nip1*njp1
        mg_off_2d_jk(1,lev) = mg_off_2d_jk(numbl,lev-1) + njp1*nkp1
        mg_off_2d_ki(1,lev) = mg_off_2d_ki(numbl,lev-1) + nkp1*nip1
        mg_off(1,lev)       = mg_off(numbl,lev-1)  + nip1*njp1*nkp1

        mg_off_1d_z(1,lev)  = mg_off_1d_z(numbl,lev-1)  + nkp1+1
        mg_off_2d_xy(1,lev) = mg_off_2d_xy(numbl,lev-1) +
     &                                                 (nip1+1)*(njp1+1)

        do 40 iblk = 1,numbl-1
          
          nip1 = mg_nijkp1(1,iblk,lev)
          njp1 = mg_nijkp1(2,iblk,lev)
          nkp1 = mg_nijkp1(3,iblk,lev)

          mg_off_1d_i(iblk+1,lev)  = mg_off_1d_i(iblk,lev)  + nip1
          mg_off_1d_j(iblk+1,lev)  = mg_off_1d_j(iblk,lev)  + njp1
          mg_off_1d_k(iblk+1,lev)  = mg_off_1d_k(iblk,lev)  + nkp1
          mg_off_2d_ij(iblk+1,lev) = mg_off_2d_ij(iblk,lev) + nip1*njp1
          mg_off_2d_jk(iblk+1,lev) = mg_off_2d_jk(iblk,lev) + njp1*nkp1
          mg_off_2d_ki(iblk+1,lev) = mg_off_2d_ki(iblk,lev) + nkp1*nip1
          mg_off(iblk+1,lev)       = mg_off(iblk,lev) +  nip1*njp1*nkp1

          mg_off_1d_z(iblk+1,lev)  = mg_off_1d_z(iblk,lev)  + nkp1+1
          mg_off_2d_xy(iblk+1,lev) = mg_off_2d_xy(iblk,lev) +
     &                                                 (nip1+1)*(njp1+1)
 40     continue
 35   continue

C-----------------------------------------------------------------------
C OFFSET FOR COARSE-TO-FINE GRID POINTS CORRESPONDENCE
C-----------------------------------------------------------------------

C i,j,k OFFSETS

      do 60 lev = 1,nummgl-1
        do 65 iblk = 1,numbl
          if (mg_nijkp1(1,iblk,lev) .eq. mg_nijkp1(1,iblk,lev+1)) then
            mg_ijk_mult(1,iblk,lev) = 1
          else
            mg_ijk_mult(1,iblk,lev) = 2
          end if
          if (mg_nijkp1(2,iblk,lev) .eq. mg_nijkp1(2,iblk,lev+1)) then
            mg_ijk_mult(2,iblk,lev) = 1
          else
            mg_ijk_mult(2,iblk,lev) = 2
          end if
          if (mg_nijkp1(3,iblk,lev) .eq. mg_nijkp1(3,iblk,lev+1)) then
            mg_ijk_mult(3,iblk,lev) = 1
          else
            mg_ijk_mult(3,iblk,lev) = 2
          end if
 65     continue
 60   continue

      return
      end
