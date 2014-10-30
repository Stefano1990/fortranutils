!-------------------------
! read in block dimensions
!-------------------------

      subroutine read_bl_size(dev)

      integer      dev

C PARAMETERS

#include "dim.inc"
#include "comm.inc"

C LOCAL

      integer      iblk, ierr
      character*32 filename

      logical available

#include "cb.comm.inc"
#include "cb.lgrid.inc"

C-----------------------------------------------------------------------

C???? TEMPORARILY ASSUME THE FILE STORING THE x-COORDINATE IS CALLED
C     x.bin

      filename = 'x.bin'

      if (master) then

        inquire(file=filename, exist=available)
        if (.not. available) then
          write(*,'(''*** Error: File '',a,'' not found'')') filename
          ierr = 2
          return
        else
c         open(unit=dev, file=filename, form='unformatted')
C???? FOR DEBUGGING, USE ASCII FILE
          open(unit=dev, file=filename)
        end if

        read(dev,*) numbl
        do iblk = 1,numbl
          read(dev,*,err=10) 
     &      bl_nijkp1(1,iblk),bl_nijkp1(2,iblk),bl_nijkp1(3,iblk)
        end do
        close(dev)

      end if

      return

 10   continue
      write(*,'(''*** Error: Reading block dimension failed'')')
      ierr = 3

      return
      end
