C=======================================================================
C FILE 3-D VARIABLES
C=======================================================================

      subroutine file_3(dev,filename,filefrmt,f,ilev,mg_or_sg,ierr)

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

C INPUT/OUTPUT

      integer       dev, ilev, mg_or_sg, filefrmt, ierr
      character*(*) filename
      real          f(*)

C LOCAL

#ifdef WRITE_DOUBLE
      real          work(maxbl2*maxijk3)
#else
      real*4        work(maxbl2*maxijk3)
#endif

      integer   iblk, ilev1, iblk0
      integer   i, j, k, ii, i0
      integer   nip1, njp1 ,nkp1
      integer*4 numbl1, itemp(3)
      integer   jo, ko, jkoff, koff, poff, itemp2(3,maxbl2)
      integer   status(MPI_STATUS_SIZE), sreq(maxbl2*3)
      integer   mpierr
      logical   isopen, available, file_bin


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C MASTER RECEIVES DATA FROM SLAVES AND FILES THE DATA
C-----------------------------------------------------------------------

      sreq = MPI_REQUEST_NULL
      if (mg_or_sg.eq.1) then
        ilev1 = ilev
      else
        ilev1 = 1
      end if

      if (master) then

C OPEN FILE

        inquire(unit=dev, opened=isopen)
        if (isopen) then
          write(*,'(''ERROR: Device number '',i3,
     &      '' already in use'')') dev
          ierr = 1
          return
        end if

        inquire(file=filename, exist=available)
        if (available) then
          write(*,'(''*** Writing file '',a)') filename
        else
          write(*,'(''*** Writing file '',a)') filename
        end if
        file_bin = filefrmt .eq. 0
        if (file_bin) then
          open(unit=dev, file=filename,form='unformatted')
        else
          open(unit=dev, file=filename)
        end if

C FILE BLOCK DIMENSION

        if (file_bin) then
          numbl1 = numbl_total
          write(dev) numbl1
        else
          write(dev,*) numbl_total
        endif

        iblk = 0
        do iblk0 = 1,numbl_total
          if(bl2pe(iblk0).eq.pe)then
            iblk = iblk + 1
            itemp(1) = mg_nijkp1(1,iblk,ilev1)
            itemp(2) = mg_nijkp1(2,iblk,ilev1)
            itemp(3) = mg_nijkp1(3,iblk,ilev1)
          else
            call MPI_RECV(itemp, 3, MPI_INTEGER4, bl2pe(iblk0),
     &      iblk0-1, MPI_COMM_WORLD, status, mpierr)
          end if
          if (file_bin) then
            write(dev) itemp(1),itemp(2),itemp(3)
          else
            write(dev,*) (itemp(i),i=1,3)
          end if
        end do


C FILE DATA

        iblk = 0
        do iblk0 = 1,numbl_total
          if(bl2pe(iblk0).eq.pe)then
            iblk = iblk + 1

            nip1 = mg_nijkp1(1,iblk,ilev1)
            njp1 = mg_nijkp1(2,iblk,ilev1)
            nkp1 = mg_nijkp1(3,iblk,ilev1)
            if (mg_or_sg.eq.1) then
              poff = mg_off(iblk,ilev1)
            else
              poff = sg_off(iblk)
            end if
            jo = nip1
            ko = nip1*njp1

            if (file_bin) then
              do k = 1,nkp1
                koff = ko*(k-1) + poff
                do j = 1,njp1
                  jkoff = jo*(j-1) + koff
                  do i = 1,nip1
                    work(i) = f(i+jkoff)
                  end do
                  write(dev,err=1000) (work(i),i=1,nip1)
                end do
              end do
            else
              do k = 1,nkp1
                koff = ko*(k-1) + poff
                do j = 1,njp1
                  jkoff = jo*(j-1) + koff
                  write(dev,'(8(1x,g15.8))',err=1000)
     &              (f(i+jkoff),i=1,nip1)
                end do
              end do
            end if

C         Data located on different pe
          else
            call MPI_RECV(itemp2, 3, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &      maxbl+iblk0-1, MPI_COMM_WORLD, status, mpierr)
            call MPI_RECV(work, itemp2(1,1)*itemp2(2,1)*itemp2(3,1),
     &      REAL_TYPE_MPI_WRITE, bl2pe(iblk0), 2*maxbl+iblk0-1,
     &      MPI_COMM_WORLD, status, mpierr)
            do k = 1,itemp2(3,1)
              do j = 1,itemp2(2,1)
                if(file_bin)then
                  write(dev,err=1000)(work((k-1)*itemp2(2,1)*itemp2(1,1)
     &            +(j-1)*itemp2(1,1)+i),i=1,itemp2(1,1))
                else
                  write(dev,'(8(1x,g15.8))',err=1000)
     &            (work((k-1)*itemp2(2,1)*itemp2(1,1)
     &            +(j-1)*itemp2(1,1)+i),i=1,itemp2(1,1))
                endif
              end do
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

            itemp(1) = mg_nijkp1(1,iblk,ilev1)
            itemp(2) = mg_nijkp1(2,iblk,ilev1)
            itemp(3) = mg_nijkp1(3,iblk,ilev1)
            call MPI_ISEND(itemp, 3, MPI_INTEGER4, 0,
     &        iblk0-1, MPI_COMM_WORLD, sreq(3*(iblk-1)+1),
     &        mpierr)
c           call MPI_SEND(itemp, 3, MPI_INTEGER4, 0,
c    &        iblk0-1, MPI_COMM_WORLD,
c    &        mpierr)
          endif
        enddo

        ii = 0
        iblk=0
        do iblk0 = 1,numbl_total
          if(pe.eq.bl2pe(iblk0))then
            iblk=iblk+1

            itemp2(1,iblk) = mg_nijkp1(1,iblk,ilev1)
            itemp2(2,iblk) = mg_nijkp1(2,iblk,ilev1)
            itemp2(3,iblk) = mg_nijkp1(3,iblk,ilev1)
            call MPI_ISEND(itemp2(1,iblk), 3, INTEGER_TYPE_MPI, 0,
     &        maxbl+iblk0-1, MPI_COMM_WORLD,
     &        sreq(3*(iblk-1)+2), mpierr)
c           call MPI_SEND(itemp2(1,iblk), 3, INTEGER_TYPE_MPI, 0,
c    &        maxbl+iblk0-1, MPI_COMM_WORLD,
c    &        mpierr)

            if (mg_or_sg.eq.1) then
              poff = mg_off(iblk,ilev1)
            else
              poff = sg_off(iblk)
            end if
            i0 = ii + 1
            jo = itemp2(1,iblk)
            ko = itemp2(1,iblk)*itemp2(2,iblk)
            do k = 1,itemp2(3,iblk)
              koff = ko*(k-1) + poff
              do j = 1,itemp2(2,iblk)
                jkoff = jo*(j-1) + koff
                do i = 1,itemp2(1,iblk)
                  ii = ii + 1
                  work(ii) = f(i+jkoff)
                end do
              end do
            end do
            call MPI_ISEND(work(i0),
     &        itemp2(1,iblk)*itemp2(2,iblk)*itemp2(3,iblk),
     &        REAL_TYPE_MPI_WRITE, 0, 2*maxbl+iblk0-1, MPI_COMM_WORLD,
     &        sreq(3*(iblk-1)+3), mpierr)
c           call MPI_SEND(work(i0),
c    &        itemp2(1,iblk)*itemp2(2,iblk)*itemp2(3,iblk),
c    &        MPI_REAL4, 0, 2*maxbl+iblk0-1, MPI_COMM_WORLD,
c    &        mpierr)

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
      write(*,'(''*** Error: Filing file '',a,'' failed'')') filename

      return
      end
