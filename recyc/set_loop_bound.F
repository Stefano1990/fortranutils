C=======================================================================
C SET LOOP INDEX BOUNDS
C=======================================================================

C NOTE: ACCORDING TO ANDY'S PLAN, ALL INDICES OF LOOPS OVER GRID POINTS
C   SHOULD START FROM 2. BUT IT IS NOT CLEAR YET TO ME HOW THIS CAN BE
C   DONE WITHOUT LARGELY MODIFYING THE PARTS ADOPTED FROM PACE, WHICH
C   IS ALWAYS CUMBERSOME. SO HERE THE OLD VERSION OF THE CODE (SINGLE
C   BLOCK) IS FOLLOWED

      subroutine set_loop_bound(ivar,ilev1,ilev2,iblk1,iblk2,
     &  choice,ijkse)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS

#include "dim.inc"
#include "les.inc"
#include "bc_type.inc"
#include "data_type.inc"
#include "comm.inc"

C INPUT

      integer ivar, ilev1, ilev2, iblk1, iblk2
      integer choice
      integer ijkse(6)

C LOCAL

      integer ilev, iblk
      integer nip1, njp1, nkp1
      integer istr, jstr, kstr
      integer ierr
      real    red
#ifdef DEBUG
      real    npts_all2(maxbl)
#endif

C INCLUDE

#include "cb.lcontrol.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do 10 ilev = ilev1,ilev2

        if(choice.eq.2)npts_all(ilev) = 0

        do 15 iblk = iblk1,iblk2

          nip1 = mg_nijkp1(1,iblk,ilev)
          njp1 = mg_nijkp1(2,iblk,ilev)
          nkp1 = mg_nijkp1(3,iblk,ilev)

C SET LOWER BOUNDS

          if (ivar.eq.SV_U.and.bl_bc(1,iblk).ne.BCT_PERIODIC .and. 
     &        bl_bc(1,iblk).ne.BCT_INTERFACE) then
            istr = 3
          else
            istr = 2
          end if
          if (ivar.eq.SV_V.and.bl_bc(3,iblk).ne.BCT_PERIODIC .and. 
     &        bl_bc(3,iblk).ne.BCT_INTERFACE) then
            jstr = 3
          else
            jstr = 2
          end if
          if (ivar.eq.SV_W.and.bl_bc(5,iblk).ne.BCT_PERIODIC .and. 
     &        bl_bc(5,iblk).ne.BCT_INTERFACE) then
            kstr = 3
          else
            kstr = 2
          end if

C SET UPPER BOUNDS

          if (choice .ge. 1) then
            ijkse(1) = istr
            ijkse(2) = nip1 - 1
            ijkse(3) = jstr
            ijkse(4) = njp1 - 1
            ijkse(5) = kstr
            ijkse(6) = nkp1 - 1
          end if
          if (choice .eq. 2) then
            nijkse(1,iblk,ilev) = istr
            nijkse(2,iblk,ilev) = nip1 - 1
            nijkse(3,iblk,ilev) = jstr
            nijkse(4,iblk,ilev) = njp1 - 1
            nijkse(5,iblk,ilev) = kstr
            nijkse(6,iblk,ilev) = nkp1 - 1
            npts_all(ilev) = npts_all(ilev) + (nip1 - istr)
     &                      * (njp1 - jstr) * (nkp1 - kstr)
#ifdef DEBUG
            npts_all2(iblk) = (nip1-istr) * (njp1-jstr) * (nkp1-kstr)
#endif
          end if


 15     continue

        if (choice .eq. 2) then
#ifdef DEBUG
          call usr_sum_double (npts_all2, npts_all(ilev), numbl_total,
     &    REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
          red = npts_all(ilev)
          call MPI_ALLREDUCE(red, npts_all(ilev), 1, REAL_TYPE_MPI,
     &    MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        end if

 10   continue

      return
      end
