C=======================================================================
C  READ 3-D QUANTITIES
C=======================================================================

      subroutine read_3(dev,filename,filefrmt,f,work,ierr)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

C INCLUDE FILES
#include "cb.lgrid.inc"
#include "cb.comm.inc"

C INPUT/OUTPUT

      integer       dev, filefrmt, ierr
      character*(*) filename
      real          f(*)
      real          work(*)
#ifdef READ_DOUBLE
      real          work2(maxijk1)
#else
      real*4        work2(maxijk1)
#endif

C LOCAL

      integer   iblk, iblk2, numbl_rd, ierr1
      integer   i, j, k
      integer   nip1, njp1 ,nkp1
      integer*4 numbl_rd1, nip11, njp11 ,nkp11
      integer   jo, ko, jkoff, koff, poff
      logical   isopen, available, read_bin
      integer   lbuf, lname
      parameter    (lbuf=32)
      character*(lbuf) buffer

      integer ipe, mpierr, status(MPI_STATUS_SIZE)
      integer ntota, nrecs, rreq(maxbl*maxijk1**2)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C MASTER DOES I/O
C-----------------------------------------------------------------------

      if (master) then

C CALCULATE THE LENGTH OF THE FILE NAME

        call extr_string(1,buffer,lname,filename,lbuf,ierr)

C OPEN FILE

        inquire(unit=dev, opened=isopen)
        if (isopen) then
          write(*,'(''*** Error: Logical device number '',i3,
     &      '' is already in use'')') dev
          ierr = 1
          return
        end if

        ierr1 = 0
        inquire(file=filename(1:lname), exist=available)
        if (available) then
          write(*,'(''Reading '',a)') filename(1:lname)
        else
          write(*,'('' '',a,'' not found.'')')
     &      filename(1:lname)
          ierr1 = -1
        end if

        if (ierr1 .eq. 0) then

          read_bin = filefrmt .eq. 0
          if (read_bin) then
            open(unit=dev, file=filename(1:lname), form='unformatted')
          else
            open(unit=dev, file=filename(1:lname))
          end if

C READ BLOCK DIMENSIONS

          if (read_bin) then
            read(dev,err=1000) numbl_rd1
            numbl_rd = numbl_rd1
          else
            read(dev,*,err=1000) numbl_rd
          end if
          if (numbl_rd.ne.numbl_total) then
            write(*,'(''*** Error: '',
     &        ''Data read inconsistent with known in file '',a)')
     &        filename(1:lname)
            ierr = 3
            return
          end if

          do iblk = 1,numbl_total
            if (read_bin) then
              read(dev,err=1000) nip11,njp11,nkp11
              nip1 = nip11
              njp1 = njp11
              nkp1 = nkp11
            else
              read(dev,*,err=1000) nip1,njp1,nkp1
            end if
            if (nip1.ne.bl_nijkp1_master(1,iblk) .or.
     &          njp1.ne.bl_nijkp1_master(2,iblk) .or.
     &          nkp1.ne.bl_nijkp1_master(3,iblk)) then
              write(*,'(''*** Error: '',
     &          ''Inconsistent block dimension for block '',i3,
     &          '' in file '',a)') iblk,filename(1:lname)
              ierr = 4
              return
            end if
          end do

        end if
        do ipe = 1,numpe-1
          call MPI_SEND(ierr1, 1, INTEGER_TYPE_MPI, ipe, 1,
     &      MPI_COMM_WORLD, mpierr)
        end do

C READ DATA BLOCK BY BLOCK

C SINGLE GRID QUANTITIES (ONLY THE PRESSURE IS SOLVED BY MULTIGRID)

        iblk2 = 0
        do iblk = 1,numbl_total
          ntota = 0

          nip1 = bl_nijkp1_master(1,iblk)
          njp1 = bl_nijkp1_master(2,iblk)
          nkp1 = bl_nijkp1_master(3,iblk)
          if (pe.eq.bl2pe(iblk)) then
            iblk2 = iblk2 + 1
            if (filename(1:1).ne.'p') then
              poff = sg_off(iblk2)
            else
              poff = mg_off(iblk2,1)
            end if
          else
            poff = 0
          end if
          jo = nip1
          ko = nip1*njp1

          if (ierr1 .eq. -1. and. pe.eq.bl2pe(iblk)) then
            do k = 1,nkp1
              koff = ko*(k-1) + poff
              do j = 1,njp1
                jkoff = jo*(j-1) + koff
                do i = 1,nip1
                  f(i+jkoff) = 0.
                end do
              end do
            end do
          elseif (ierr1.ne.-1)then
            do k = 1,nkp1
              koff = ko*(k-1) + poff
              do j = 1,njp1
                jkoff = jo*(j-1) + koff
                if (read_bin) then
                  read(dev,err=1000)(work2(i),i=1,nip1)
                  do i = 1,nip1
                    work(i)=work2(i)
                  end do
                else
                  read(dev,*,err=1000)(work(i),i=1,nip1)
                end if
                if (pe.eq.bl2pe(iblk)) then
                  do i = 1,nip1
                    f(i+jkoff) = work(i)
                  end do
                else
                  call MPI_SEND(work, nip1, REAL_TYPE_MPI,
     &            bl2pe(iblk), iblk, MPI_COMM_WORLD, mpierr)
                end if 
              end do
            end do
          end if

        end do

        if (ierr1 .eq. 0) then
          close(dev)
        else if(ierr1 .lt. 0) then
          ierr1 = 0
        end if

C-----------------------------------------------------------------------
C SLAVES RECEIVE DATA FROM MASTER
C-----------------------------------------------------------------------

      else

        call MPI_RECV(ierr1, 1, INTEGER_TYPE_MPI, 0, 1, MPI_COMM_WORLD,
     &    status, mpierr)

        nrecs = 0
        iblk2 = 0

        do iblk = 1,numbl_total

          if(pe.eq.bl2pe(iblk))then

            iblk2 = iblk2 + 1

            ntota = 0
            nip1 = bl_nijkp1(1,iblk2)
            njp1 = bl_nijkp1(2,iblk2)
            nkp1 = bl_nijkp1(3,iblk2)
            if (filename(1:1).ne.'p') then
              poff = sg_off(iblk2)
            else
              poff = mg_off(iblk2,1)
            end if
            jo = nip1
            ko = nip1*njp1

            if (ierr1 .eq. -1) then
              do k = 1,nkp1
                koff = ko*(k-1) + poff
                do j = 1,njp1
                  jkoff = jo*(j-1) + koff
                    do i = 1,nip1
                    f(i+jkoff) = 0.
                  end do
                end do
              end do
            else
              do k = 1,nkp1
                koff = ko*(k-1) + poff
                do j = 1,njp1
                  jkoff = jo*(j-1) + koff
                  nrecs = nrecs + 1
                  call MPI_IRECV(f(jkoff+1), nip1, REAL_TYPE_MPI,
     &            0, iblk, MPI_COMM_WORLD, rreq(nrecs), mpierr)
                end do
              end do
            end if

          end if

        end do

        do i = 1,nrecs
          call MPI_WAIT(rreq(i), status, ierr)
        end do

      end if

      return

 1000 ierr = 5
      write(*,'(''*** Error: Reading file '',a,'' failed'')')
     &  filename(1:lname)

      return
      end
