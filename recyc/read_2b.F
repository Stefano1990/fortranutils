C=======================================================================
C FILE 2-D VARIABLES
C=======================================================================

C NOTE: THE 2-D QUANTITY IS IN A PLANE DETERMINED BY d2-d3. THE ORDERING
C OF THE QUANTITY IS IN NATURAL ORDER WITH d2 COUNTED FIRST

      subroutine read_2b(dev,filename,d2,d3,f,mg_or_sg,rf,read_status,
     &  ilev,work,ierr)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

C INCLUDE

#include "cb.lgrid.inc"
#include "cb.comm.inc"
#include "cb.lcontrol.inc"

C INPUT/OUTPUT

      integer       d2, d3
      integer       rf, read_status
      integer       dev, mg_or_sg, ilev, ierr
      character*(*) filename
      real          f(*), work(*)

C LOCAL

      integer   iblk, iblk2
      integer   d1
      integer   i, j, k
      integer   njp1 ,nkp1
      integer*4 numbl1, nip11, njp11 ,nkp11
      integer   koff, poff
      logical   isopen, available, binary_read

      integer ipe, mpierr, status(MPI_STATUS_SIZE)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      d1 = 6 - (d2+d3)
      binary_read = rf .eq. 0

C-----------------------------------------------------------------------
C MASTER READ DATA
C-----------------------------------------------------------------------

      if (master) then

        if (read_status .eq. 0) then

C OPEN FILE

          inquire(unit=dev, opened=isopen)
          if (isopen) then
            write(*,'(''*** Error: Logical device number '',i3,
     &        '' is already in use'')') dev
            ierr = 1
          end if

          inquire(file=filename, exist=available)
          if (available) then
            write(*,'(''*** Info: Reading file '',a)') filename
          else
            write(*,'(''*** Info: File '',a,
     &        '' not found; set to zero'')') filename
            ierr = -1
          end if
          if (ierr .gt. 0) return

          if (ierr .eq. 0) then
            if (binary_read) then
              open(unit=dev, file=filename, form='unformatted')
            else
              open(unit=dev, file=filename)
            end if
          end if

C READ BLOCK DIMENSION

          if (ierr .eq. 0) then
            if (binary_read) then
              read(dev) numbl1
              do iblk = 1,numbl1
                read(dev) nip11, njp11, nkp11
              end do
            else
              read(dev,*) numbl1
              do iblk = 1,numbl1
                read(dev,*) nip11, njp11, nkp11
              end do
            end if
          end if

        end if

C TRANSMIT IERR
        do ipe = 2,numpe
          call MPI_SEND(ierr, 1, INTEGER_TYPE_MPI, ipe-1, 1,
     &      MPI_COMM_WORLD, mpierr)
        end do


C READ DATA

        iblk2 = 0
        do iblk = 1,numbl_total

c         njp1 = mg_nijkp1(d2,iblk,ilev)
c         nkp1 = mg_nijkp1(d3,iblk,ilev)
C         Some communication would be necessary here for ilev =/= 1.
C         Ignoring this case as ilev is always 1. But note that some
C         changes would be required if it wasn't.
          njp1 = bl_nijkp1_master(d2,iblk)
          nkp1 = bl_nijkp1_master(d3,iblk)

          if (pe.eq.bl2pe(iblk)) then
            iblk2 = iblk2 + 1
            if (mg_or_sg .eq. 1) then
              if (d1 .eq. 1) then
                poff = mg_off_2d_jk(iblk2,ilev)
              else if (d1.eq.2) then
                poff = mg_off_2d_ki(iblk2,ilev)
              else
                poff = mg_off_2d_ij(iblk2,ilev)
              end if
            else
              if (d1 .eq. 1) then
                poff = sg_off_2d_jk(iblk2)
              else if (d1.eq.2) then
                poff = sg_off_2d_ki(iblk2)
              else
                poff = sg_off_2d_ij(iblk2)
              end if
            end if
          end if

          if (ierr .eq. 0) then
            do k = 1,nkp1
              if (binary_read) then
                read(dev,err=1000) (work(j),j=1,njp1)
              else
                read(dev,'(8(1x,g15.8))',err=1000) (work(j),j=1,njp1)
              end if
              if (pe.eq.bl2pe(iblk)) then
                koff = njp1*(k-1) + poff
                do j = 1,njp1
                  f(j+koff) = work(j)
                end do
              else
                call MPI_SEND(work, njp1, REAL_TYPE_MPI,
     &          bl2pe(iblk), 2, MPI_COMM_WORLD, mpierr)
              end if 
            end do
          else if (ierr .lt. 0 .and. pe.eq.bl2pe(iblk)) then
            do k = 1,nkp1
              koff = njp1*(k-1) + poff
              do j = 1,njp1
                f(j+koff) = 0.
              end do
            end do
          end if

        end do

        if (read_status.eq.2 .and. ierr.eq.0) then
          close(dev)
          ierr = 0
        end if

C-----------------------------------------------------------------------
C SLAVES RECEIVE DATA FROM THE MASTER
C-----------------------------------------------------------------------

      else

        call MPI_RECV(ierr, 1, INTEGER_TYPE_MPI, 0, 1, MPI_COMM_WORLD,
     &    status, mpierr)

        do iblk = 1,numbl

          njp1 = bl_nijkp1(d2,iblk)
          nkp1 = bl_nijkp1(d3,iblk)

          if (mg_or_sg .eq. 1) then
            if (d1 .eq. 1) then
              poff = mg_off_2d_jk(iblk,ilev)
            else if (d1.eq.2) then
              poff = mg_off_2d_ki(iblk,ilev)
            else
              poff = mg_off_2d_ij(iblk,ilev)
            end if
          else
            if (d1 .eq. 1) then
              poff = sg_off_2d_jk(iblk)
            else if (d1.eq.2) then
              poff = sg_off_2d_ki(iblk)
            else
              poff = sg_off_2d_ij(iblk)
            end if
          end if

          if (ierr .eq. 0) then
            do k = 1,nkp1
              koff = njp1*(k-1) + poff
              call MPI_RECV(f(koff+1), njp1, REAL_TYPE_MPI,
     &        0, 2, MPI_COMM_WORLD, status, mpierr)
            end do
          else if (ierr .lt. 0) then
            do k = 1,nkp1
              koff = njp1*(k-1) + poff
              do j = 1,njp1
                f(j+koff) = 0.
              end do
            end do
          end if

        end do

      end if

      return

 1000 ierr = 5
      write(*,'(''*** Error: Reading file '',a,'' failed'')') filename

      return
      end

