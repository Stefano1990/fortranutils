!-------------------------------
! read in grid co-ordinates
! x and y are 2-d arrays
! z is a 1-d array
!-------------------------------

      subroutine read_xyz(dev,filename,xyz,work,ierr)

      implicit none


C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

C INPUT/OUTPUT

      integer       dev, ierr
      real          xyz(*), work(*)
      character*(*) filename

C LOCAL

      integer i, j, k, joff, poff
      integer nip1, njp1, nkp1, nip2
      integer iblk, iblk2, numbl_rd
      logical isopen, available
      integer status(MPI_STATUS_SIZE), icount


C INCLUDE

#include "cb.lgrid.inc"
#include "cb.comm.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      if (master) then

C OPEN FILE

        inquire(unit=dev, opened=isopen)
        if (isopen) then
          write(*,'(''*** Error: Logical device number '',i3,
     &      '' is already in use'')') dev
          ierr = 1
          return
        end if

        inquire(file=filename, exist=available)
        if (available) then
          write(*,'(''Reading '',a)') filename
        else
          write(*,'(''File '',a,'' not found'')') filename
          ierr = 2
          return
        end if
        open(unit=dev, file=filename)

C READ BLOCK DIMENSIONS

        read(dev,*,err=1000) numbl_rd
        if (numbl_rd.ne.numbl_total) then
          write(*,'(''*** Error: '',
     &      ''Data read inconsistent with known in file '',a)') filename
          ierr = 3
          return
        end if

        do iblk = 1,numbl_total
          read(dev,*,err=1000) nip1,njp1,nkp1
          if (nip1.ne.bl_nijkp1_master(1,iblk) .or.
     &        njp1.ne.bl_nijkp1_master(2,iblk) .or.
     &        nkp1.ne.bl_nijkp1_master(3,iblk)) then
            write(*,'(''Error: '',
     &        ''Inconsistent block dimension for block '',i3,
     &        '' in file'',a)') iblk,filename
            ierr = 4
            return
          end if
        end do

C READ COORDINATES (DUE TO THE PACE'S SPECIAL DATA SUPPLY OF THE
C COORDINATES, READ HAS TO BE DONE SEPARATELY FOR x- (y-) AND z.

        if (filename(1:1).eq.'z') then

C - READ z-COORDINATE

          iblk2 = 0
          do iblk = 1,numbl_total

            nkp1 = bl_nijkp1_master(3,iblk)
            if(pe.eq.bl2pe(iblk))then
              iblk2 = iblk2 + 1
              poff = mg_off_1d_z(iblk2,1)
              read(dev,*,err=1000) (xyz(poff+k), k=2,nkp1)
            else
              read(dev,*,err=1000) (work(k), k=1,nkp1-1)
              call MPI_SEND(work, nkp1-1, REAL_TYPE_MPI,
     &          bl2pe(iblk), 1, MPI_COMM_WORLD, ierr)
            endif
          end do

        else

C - READ x OR y-COORDINATE

          iblk2 = 0
          do iblk = 1,numbl_total
            nip1 = bl_nijkp1_master(1,iblk)
            njp1 = bl_nijkp1_master(2,iblk)
            nkp1 = bl_nijkp1_master(3,iblk)

            if(pe.eq.bl2pe(iblk))then
              iblk2 = iblk2 + 1
              poff = mg_off_2d_xy(iblk2,1)
              nip2 = nip1 + 1

              do j = 2,njp1
                joff = poff + nip2*(j-1)
                read(dev,*,err=1000) (xyz(joff+i), i=2,nip1)
              end do

            else
              do j = 2,njp1
                read(dev,*,err=1000) (work(i), i=1,nip1-1)
                call MPI_SEND(work, nip1-1, REAL_TYPE_MPI,
     &            bl2pe(iblk), j, MPI_COMM_WORLD, ierr)
              end do
            endif

          end do

        end if

C CLOSE FILE
        close(dev)

      else

C - RECEIVE z-COORDINATE
        if (filename(1:1).eq.'z') then
          do iblk = 1,numbl
            poff = mg_off_1d_z(iblk,1)
            nkp1 = mg_nijkp1(3,iblk,1)
            call MPI_RECV(xyz(poff+2), nkp1-1, REAL_TYPE_MPI,
     &        0, 1, MPI_COMM_WORLD, status, ierr)
          end do

C - RECEIVE x OR y-COORDINATE
        else
          do iblk = 1,numbl
            poff = mg_off_2d_xy(iblk,1)
            nip1 = mg_nijkp1(1,iblk,1)
            njp1 = mg_nijkp1(2,iblk,1)
            nip2 = nip1 + 1
            do j = 2,njp1
              joff = poff + nip2*(j-1)
              call MPI_RECV(xyz(joff+2), nip1-1, REAL_TYPE_MPI,
     &          0, j, MPI_COMM_WORLD, status, ierr)
            end do
          end do
        end if

      end if


      return

 1000 ierr = 5
      write(*,'(''*** Error: Reading file '',a,'' failed'')') filename

      return
      end
