C=======================================================================
C READ BLOCK TO BLOCK AND BOUNDARY INFORMATION
C=======================================================================

C ????? IT IS NOT KNOWN IN WHAT FORMAT ANDY WOULD LIKE TO SUPPLY THIS
C INFORMATION. SO EVERYTHING DONE HERE HAS TO BE BASED ON ASSUMPTIONS

      subroutine read_bl_info(ierr)

C MPI
#include "mpif.h"

C PARAMETERS

#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"
        integer      ierr


C LOCAL

      integer      iblk, iside, dev, lname
      character*32 fname, filename

      logical      available

      integer      numbl_pe(maxbl), ipe, status(MPI_STATUS_SIZE)
      integer      mpierr

C INCLUDE

#include "cb.lgrid.inc"
#include "cb.comm.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      ierr = 0

      if (master) then

        filename = 'block.dat'
        call extr_string(1,fname,lname,filename,32,ierr)

        do dev = 10,99
          inquire(unit=dev, opened=available)
          if (.not.available) goto 10
        end do
 10     continue

        inquire(file=filename, exist=available)
        if (available) then
          write(*,'(''Reading '',a)') fname(1:lname)
          open(unit=dev, file=filename)
        else
          write(*,'(''ERROR:  File '',a,'' is missing'')')
     &      fname(1:lname)
          ierr = 1
          return
        end if

C INITIALIZE ALL INFORMATION ARRAYS

        do iblk = 1,maxbl
          bl2pe(iblk) = 0
          numbl_pe(iblk) = 0
          do iside = 1,6
            bl2bl(iside,iblk)    = 0
            bl2side(iside,iblk)  = 0
            bl_bc(iside,iblk)    = 0
            bl_split(iside,iblk) = 0
          end do
        end do

C BLOCK ASSIGMENT, BLOCK TO BLOCK CONNECTION, AND BLOCK SIDE CONNECTION
C INFORMATION

        read(dev,*)
        do iblk = 1,numbl
          read(dev,*,err=1000)
     &     bl2pe(iblk),(bl2bl(iside,iblk),bl2side(iside,iblk),iside=1,6)
          numbl_pe(bl2pe(iblk)+1)=numbl_pe(bl2pe(iblk)+1)+1
        end do

C BOUNDARY CONDITION

        read(dev,*)
        do iblk = 1,numbl
          read(dev,*,err=1005) (bl_bc(iside,iblk),iside=1,6)
        end do

        do ipe=1,numpe-1
          CALL MPI_SEND(bl2pe, maxbl, INTEGER_TYPE_MPI,
     &    ipe, 1, MPI_COMM_WORLD, mpierr)
          CALL MPI_SEND(bl2bl, 6*maxbl, INTEGER_TYPE_MPI,
     &    ipe, 2, MPI_COMM_WORLD, mpierr)
          CALL MPI_SEND(bl2side, 6*maxbl, INTEGER_TYPE_MPI,
     &    ipe, 3, MPI_COMM_WORLD, mpierr)
          CALL MPI_SEND(numbl_pe(ipe+1), 1, INTEGER_TYPE_MPI,
     &    ipe, 4, MPI_COMM_WORLD, mpierr)
        end do

        do iblk = 1,numbl
          if(pe.ne.bl2pe(iblk))then
            CALL MPI_SEND(bl_bc(1,iblk), 6, INTEGER_TYPE_MPI,
     &      bl2pe(iblk), 5, MPI_COMM_WORLD, mpierr)
          end if
        end do

C OUTFLOW MASS

!        read(dev,*)
!        do iblk = 1,numbl
!          read(dev,*,err=1010) (bl_flout(iside,iblk),iside=1,6)
!        end do

!        do iblk = 1,numbl
!          if(pe.ne.bl2pe(iblk))then
!            CALL MPI_SEND(bl_flout(1,iblk), 6, REAL_TYPE_MPI,
!     &      bl2pe(iblk), 6, MPI_COMM_WORLD, mpierr)
!          end if
!        end do

C OUTFLOW SPLIT

        read(dev,*)
        do iblk = 1,numbl
          read(dev,*,err=1015) (bl_split(iside,iblk),iside=1,6)
        end do

        do iblk = 1,numbl
          if(pe.ne.bl2pe(iblk))then
            CALL MPI_SEND(bl_split(1,iblk), 6, REAL_TYPE_MPI,
     &      bl2pe(iblk), 7, MPI_COMM_WORLD, mpierr)
          end if
        end do



C PRESSURE JUMP

!        read(dev,*)
!        do iblk = 1,numbl
!          read(dev,*,err=1020) (bl_pjump(iside,iblk),iside=1,6)
!        end do

!        do iblk = 1,numbl
!          if(pe.ne.bl2pe(iblk))then
!            CALL MPI_SEND(bl_pjump(1,iblk), 6, REAL_TYPE_MPI,
!     &      bl2pe(iblk), 8, MPI_COMM_WORLD, mpierr)
!          end if
!        end do

C       Set own blocks
        iitemp = 1
        do iblk = 1,numbl
          if(pe.eq.bl2pe(iblk)) then
            do iside = 1,6
              bl_bc(iside,iitemp) = bl_bc(iside,iblk)
!              bl_flout(iside,iitemp) = bl_flout(iside,iblk)
              bl_split(iside,iitemp) = bl_split(iside,iblk)
!              bl_pjump(iside,iitemp) = bl_pjump(iside,iblk)
            end do
            iitemp = iitemp + 1
          endif
        end do

        close(dev)

C       Henceforth forget about blocks that do not belong to us
C       unless we need to distribute or collect info to or from
C       other processes. numbl now refers only to the local
C       blocks; numbl_total is the total number of blocks over
C       all processes.

        numbl_total=numbl
        numbl=numbl_pe(pe+1)

        do ipe=1,numpe-1
          CALL MPI_SEND(numbl_total, 1, INTEGER_TYPE_MPI,
     &    ipe, 9, MPI_COMM_WORLD, mpierr)
        end do

        return

      else

        CALL MPI_RECV(bl2pe, maxbl, INTEGER_TYPE_MPI, 0, 1,
     &  MPI_COMM_WORLD, status, mpierr)

        CALL MPI_RECV(bl2bl, 6*maxbl, INTEGER_TYPE_MPI, 0, 2,
     &  MPI_COMM_WORLD, status, mpierr)

        CALL MPI_RECV(bl2side, 6*maxbl, INTEGER_TYPE_MPI, 0, 3,
     &  MPI_COMM_WORLD, status, mpierr)

        CALL MPI_RECV(numbl, 1, INTEGER_TYPE_MPI, 0, 4,
     &  MPI_COMM_WORLD, status, mpierr)

        do iblk = 1,numbl
          CALL MPI_RECV(bl_bc(1,iblk), 6, INTEGER_TYPE_MPI, 0, 5,
     &    MPI_COMM_WORLD, status, mpierr)
        end do

!        do iblk = 1,numbl
!          CALL MPI_RECV(bl_flout(1,iblk), 6, REAL_TYPE_MPI,
!     &    0, 6, MPI_COMM_WORLD, status, mpierr)
!        end do

        do iblk = 1,numbl
          CALL MPI_RECV(bl_split(1,iblk), 6, REAL_TYPE_MPI,
     &    0, 7, MPI_COMM_WORLD, status, mpierr)
        end do

!        do iblk = 1,numbl
!          CALL MPI_RECV(bl_pjump(1,iblk), 6, REAL_TYPE_MPI,
!     &    0, 8, MPI_COMM_WORLD, status, mpierr)
!        end do

        CALL MPI_RECV(numbl_total, 1, INTEGER_TYPE_MPI, 0, 9,
     &  MPI_COMM_WORLD, status, mpierr)

        return

      end if


 1000 write(*,'(''ERROR: '',
     &          ''Block connectivity read in failed'')')
      ierr = 2
      return

 1005 write(*,'(''ERROR: '',
     &          ''Block boundary condition read in failed'')')
      ierr = 3
      return

 1010 write(*,'(''ERROR: '',
     &          ''Reading block mass flux condition failed'')')
      ierr = 4
      return

 1015 write(*,'(''*** Error: '',
     &          ''Reading block mass flux split condition failed'')')
      ierr = 5
      return

 1020 write(*,'(''*** Error: '',
     &          ''Reading block pressure jump condition failed'')')
      ierr = 6
      return

      end

