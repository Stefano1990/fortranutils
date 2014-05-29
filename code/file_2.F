C=======================================================================
C FILE 2-D VARIABLES
C=======================================================================

      subroutine file_2(dev,filename,d1,f,mg_or_sg,ilev,ierr)

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

      integer       dev, d1, mg_or_sg, ilev, ierr
      character*(*) filename
      real          f(*)

C LOCAL

      integer iblk, iblk0, itemp(3), itemp2(2,maxbl2)
      integer i, j, k, ii, i0
      integer njp1 ,nkp1
      integer koff, poff
      integer d2, d3
      integer status(MPI_STATUS_SIZE), sreq(3*maxbl2)
      integer mpierr
      logical isopen, available
      real    rtemp(maxbl2*maxijk1**2)


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C MASTER RECEIVES DATA FROM SLAVES AND FILES THE DATA
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
          write(*,'(''*** Info: Overwriting file '',a)') filename
        else
          write(*,'(''*** Info: Creating file '',a)') filename
        end if
        open(unit=dev, file=filename)

C FILE BLOCK DIMENSION

        write(dev,*) numbl_total
        iblk=0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk=iblk+1
            write(dev,*) (bl_nijkp1(i,iblk),i=1,3)
          else
            call MPI_RECV(itemp, 3, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &      iblk0-1, MPI_COMM_WORLD, status, mpierr)
            write(dev,*) (itemp(i),i=1,3)
          endif
        end do

C FILE DATA

        d2 = cyclic(1,d1)
        d3 = cyclic(2,d1)

        iblk=0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk=iblk+1

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

            do k = 1,nkp1
              koff = njp1*(k-1) + poff
              write(dev,'(8(1x,g15.8))',err=1000) (f(j+koff),j=1,njp1)
            end do

C         Data located on different pe
          else
            call MPI_RECV(itemp, 2, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &      maxbl+iblk0-1, MPI_COMM_WORLD, status, mpierr)
            call MPI_RECV(rtemp, itemp(1)*itemp(2), REAL_TYPE_MPI,
     &      bl2pe(iblk0), 2*maxbl+iblk0-1, MPI_COMM_WORLD, status,
     &      mpierr)
            do k = 1,itemp(2)
              write(dev,*) (rtemp((k-1)*itemp(1)+j),j=1,itemp(1))
            end do

          end if

        end do

        close(dev)

C-----------------------------------------------------------------------
C SLAVES SEND DATA TO THE MASTER
C-----------------------------------------------------------------------

      else


        iblk=0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk=iblk+1
            itemp(1) = bl_nijkp1(1,iblk)
            itemp(2) = bl_nijkp1(2,iblk)
            itemp(3) = bl_nijkp1(3,iblk)
            call MPI_ISEND(itemp, 3, INTEGER_TYPE_MPI, 0,
     &        iblk0-1, MPI_COMM_WORLD, sreq(3*(iblk-1)+1),
     &        mpierr)
          endif
        enddo

        ii = 0
        i0 = 0
        iblk=0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk=iblk+1

            d2 = cyclic(1,d1)
            d3 = cyclic(2,d1)
            itemp2(1,iblk) = mg_nijkp1(d2,iblk,ilev)
            itemp2(2,iblk) = mg_nijkp1(d2,iblk,ilev)
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
     &        REAL_TYPE_MPI, 0, 2*maxbl+iblk0-1, MPI_COMM_WORLD,
     &        sreq(3*(iblk-1)+3), mpierr)

          endif
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
      write(*,'(''*** Error: Filing file '',a,'' failed'')') filename

      return
      end
