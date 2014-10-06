!----------------------------------------------
! routines to open and close files used for
! transient flow outputs
!----------------------------------------------


      subroutine open_files
  
  
      implicit none
#include "dim.inc"
#include "const.inc"
#include "les.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "comm.inc"
#include "cb.comm.inc"

       integer iblk, dm,fm, gm, hm, mm, sm

      do iblk = 1, numbl
        dm = 1000 + iblkl2g(iblk)
        fm = 3000 + iblkl2g(iblk)
        gm = 4000 + iblkl2g(iblk)
        mm = 5000 + iblkl2g(iblk)
        sm = 6000 + iblkl2g(iblk)
!        write(*,*) pe, iblk, iblkl2g(iblk),dm
      enddo
 
       open(unit=dm)
       open(unit=fm)
       open(unit=gm)
!       open(unit=hm)
       open(unit=mm)
       open(unit=sm)


!       if(master) then
!       open(unit=1982,file='mean.dat')
!       open(unit=1983,file='stress.dat')
!       endif
   
      write(dm,*) 'VARIABLES = "x" "y" "u" "v" "w" "xi" "drho" "evisc"'
      write(fm,*) 'VARIABLES = "x" "y" "u" "v" "w" "xi" "drho" "d2rho"'
      write(gm,*) 'VARIABLES = "x" "y" "wx" "wz"'
!      write(hm,*) 'VARIABLES = "x" "z" "drho" "d2rho"'



      return
      end
  
!---------------------------------------------
