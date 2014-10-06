C=======================================================================
C FILE WALL SHEAR STRESSES
C=======================================================================

      subroutine read_twall(tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &                      tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &                      tw_wfk,tw_efk,tw_sfk,tw_nfk,
     &                      work,ierr)

        implicit none

#include "dim.inc"

#include "cb.control.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "cb.lgrid.inc"

      integer ierr
      real    tw_sfi(*), tw_nfi(*), tw_lfi(*), tw_rfi(*),
     &        tw_wfj(*), tw_efj(*), tw_lfj(*), tw_rfj(*),
     &        tw_wfk(*), tw_efk(*), tw_sfk(*), tw_nfk(*),
     &        work

      integer binary_read, ierr1

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      ierr1 = 0
#ifdef READ_DOUBLE
      binary_read = 0
#else
      binary_read = 1
#endif

      call read_2b(ndn,'wallstress',3,1,tw_sfi,2,binary_read,0,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',3,1,tw_nfi,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',1,2,tw_lfi,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',1,2,tw_rfi,2,binary_read,1,1,work,
     & ierr1)

      call read_2b(ndn,'wallstress',3,2,tw_wfj,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',3,2,tw_efj,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',1,2,tw_lfj,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',1,2,tw_rfj,2,binary_read,1,1,work,
     & ierr1)

      call read_2b(ndn,'wallstress',3,2,tw_wfk,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',3,2,tw_efk,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',3,1,tw_sfk,2,binary_read,1,1,work,
     & ierr1)
      call read_2b(ndn,'wallstress',3,1,tw_nfk,2,binary_read,2,1,work,
     & ierr1)

      if (ierr1 .gt. 0) ierr = ierr1

      return
      end
