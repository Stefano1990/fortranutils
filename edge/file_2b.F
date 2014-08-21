C=======================================================================
C FILE 2-D VARIABLES (QUANTITIES)
C=======================================================================

C NOTE: THE 2-D QUANTITY IS IN A PLANE DETERMINED BY d2-d3. THE ORDERING
C OF THE QUANTITY IS IN NATURAL ORDER WITH d2 COUNTED FIRST

      subroutine file_2b(dev,filename,d2,d3,f,mg_or_sg,wf,write_status,
     &  ilev,ierr)

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
      integer       wf, write_status
      integer       dev, mg_or_sg, ilev, ierr
      character*(*) filename
      real          f(*)

C LOCAL

      integer   iblk, iblk0
      integer   d1
      integer   i, j, k, ii, i0
      integer   njp1 ,nkp1
      integer*4 numbl1, itemp(3,maxbl2)
      integer   koff, poff, itemp2(2,maxbl2)
      integer   status(MPI_STATUS_SIZE), sreq(3*maxbl2)
      integer   mpierr
      logical   isopen, available, binary_write
#ifdef WRITE_DOUBLE
      real      rtemp(maxbl2*maxijk1**2)
#else
      real*4    rtemp(maxbl2*maxijk1**2)
#endif

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      d1 = 6 - (d2+d3)
c     d1 = 6 - d2+d3
      binary_write = wf .eq. 0

C-----------------------------------------------------------------------
C MASTER RECEIVES DATA FROM SLAVES AND FILES THE DATA
C-----------------------------------------------------------------------

      if (master) then

        if (write_status .eq. 0) then

C OPEN FILE

          inquire(unit=dev, opened=isopen)
          if (isopen) then
            write(*,'(''*** Error: Logical device number '',i3,
     &        '' is already in use'')') dev
            ierr = 1
            return
          end if

          inquire(file=filename, exist=available)
          if (available) then
            write(*,'(''Writing'',a)') filename
          else
            write(*,'(''Writing'',a)') filename
          end if

          if (binary_write) then
            open(unit=dev, file=filename, form='unformatted')
          else
            open(unit=dev, file=filename)
          end if

C FILE BLOCK DIMENSION

          if (binary_write) then
            numbl1 = numbl_total
            write(dev) numbl1
          else
            write(dev,*) numbl_total
          endif

          iblk = 0
          do iblk0 = 1,numbl_total
            if(bl2pe(iblk0).eq.pe)then
              iblk = iblk + 1
              itemp(1,1) = bl_nijkp1(1,iblk)
              itemp(2,1) = bl_nijkp1(2,iblk)
              itemp(3,1) = bl_nijkp1(3,iblk)
            else
              call MPI_RECV(itemp(1,1), 3, MPI_INTEGER4, bl2pe(iblk0),
     &        iblk0-1, MPI_COMM_WORLD, status, mpierr)
            end if
            if (binary_write) then
              write(dev) itemp(1,1),itemp(2,1),itemp(3,1)
            else
              write(dev,*) (itemp(i,1),i=1,3)
            end if
          end do

        end if

C FILE DATA

        iblk = 0
        do iblk0 = 1,numbl_total

          if(bl2pe(iblk0).eq.pe)then
            iblk = iblk + 1

            njp1 = mg_nijkp1(d2,iblk,ilev)
            nkp1 = mg_nijkp1(d3,iblk,ilev)
            if (mg_or_sg .eq. 1) then
              if (d1.eq.1) then
                poff = mg_off_2d_jk(iblk,ilev)
              else if (d1.eq.2) then
                poff = mg_off_2d_ki(iblk,ilev)
              else
                poff = mg_off_2d_ij(iblk,ilev)
              end if
            else
              if (d1.eq.1) then
                poff = sg_off_2d_jk(iblk)
              else if (d1.eq.2) then
                poff = sg_off_2d_ki(iblk)
              else
                poff = sg_off_2d_ij(iblk)
              end if
            end if

            if (binary_write) then
              do k = 1,nkp1
                koff = njp1*(k-1) + poff
                write(dev,err=1000) (f(j+koff),j=1,njp1)
              end do
            else
              do k = 1,nkp1
                koff = njp1*(k-1) + poff
                write(dev,'(8(1x,g15.8))',err=1000) (f(j+koff),j=1,njp1)
              end do
            end if

C         Data located on different pe
          else
            call MPI_RECV(itemp2(1,1), 2, INTEGER_TYPE_MPI,
     &      bl2pe(iblk0), maxbl+iblk0-1, MPI_COMM_WORLD, status,
     &      mpierr)
            call MPI_RECV(rtemp, itemp2(1,1)*itemp2(2,1),
     &      REAL_TYPE_MPI_WRITE, bl2pe(iblk0), 2*maxbl+iblk0-1,
     &      MPI_COMM_WORLD, status, mpierr)
            do k = 1,itemp2(2,1)
              if(binary_write)then
                write(dev,err=1000) (rtemp((k-1)*itemp2(1,1)+j),
     &          j=1,itemp2(1,1))
              else
                write(dev,'(8(1x,g15.8))',err=1000)
     &          (rtemp((k-1)*itemp2(1,1)+j),j=1,itemp2(1,1))
              endif
            end do

          end if

        end do

        if (write_status .eq. 2) close(dev)

C-----------------------------------------------------------------------
C SLAVES SEND DATA TO THE MASTER
C-----------------------------------------------------------------------

      else

        if (write_status .eq. 0) then
          iblk=0
          do iblk0 = 1,numbl_total
            if(pe.eq.bl2pe(iblk0))then
              iblk=iblk+1
              itemp(1,iblk) = bl_nijkp1(1,iblk)
              itemp(2,iblk) = bl_nijkp1(2,iblk)
              itemp(3,iblk) = bl_nijkp1(3,iblk)
              call MPI_ISEND(itemp(1,iblk), 3, MPI_INTEGER4, 0,
     &          iblk0-1, MPI_COMM_WORLD, sreq(3*(iblk-1)+1),
     &          mpierr)
            endif
          enddo
        endif

        ii = 0
        i0 = 0
        iblk=0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk=iblk+1

            itemp2(1,iblk) = mg_nijkp1(d2,iblk,ilev)
            itemp2(2,iblk) = mg_nijkp1(d3,iblk,ilev)
            call MPI_ISEND(itemp2(1,iblk), 2, INTEGER_TYPE_MPI, 0,
     &        maxbl+iblk0-1, MPI_COMM_WORLD,
     &        sreq(3*(iblk-1)+2), mpierr)

            if (mg_or_sg .eq. 1) then
              if (d1.eq.1) then
                poff = mg_off_2d_jk(iblk,ilev)
              else if (d1.eq.2) then
                poff = mg_off_2d_ki(iblk,ilev)
              else
                poff = mg_off_2d_ij(iblk,ilev)
              end if
            else
              if (d1.eq.1) then
                poff = sg_off_2d_jk(iblk)
              else if (d1.eq.2) then
                poff = sg_off_2d_ki(iblk)
              else
                poff = sg_off_2d_ij(iblk)
              end if
            end if

            i0 = ii + 1
            do k = 1,itemp2(2,iblk)
              koff = itemp2(1,iblk)*(k-1) + poff
              do j = 1, itemp2(1,iblk)
                ii = ii + 1
                rtemp(ii) = f(j+koff)
              end do
            end do
            call MPI_ISEND(rtemp(i0), itemp2(1,iblk)*itemp2(2,iblk),
     &        MPI_REAL4, 0, 2*maxbl+iblk0-1, MPI_COMM_WORLD,
     &        sreq(3*(iblk-1)+3), mpierr)

          end if
        end do

        iblk = 0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk = iblk + 1
            call MPI_WAIT(sreq(3*(iblk-1)+1), status, mpierr)
            call MPI_WAIT(sreq(3*(iblk-1)+2), status, mpierr)
            call MPI_WAIT(sreq(3*(iblk-1)+3), status, mpierr)
          end if
        end do


      end if

      return

 1000 ierr = 5
      write(*,'(''Writing '',a,'' failed'')') filename

      return
      end
