C=======================================================================
C FILE WALL SHEAR STRESSES
C=======================================================================

      subroutine file_twall(tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &                      tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &                      tw_wfk,tw_efk,tw_sfk,tw_nfk)

#include "dim.inc"
#include "comm.inc"

#include "cb.control.inc"
#include "cb.comm.inc"

      real    tw_sfi(*), tw_nfi(*), tw_lfi(*), tw_rfi(*),
     &        tw_wfj(*), tw_efj(*), tw_lfj(*), tw_rfj(*),
     &        tw_wfk(*), tw_efk(*), tw_sfk(*), tw_nfk(*)

      integer binary_write, ierr

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      ierr = 0
#ifdef WRITE_DOUBLE
      binary_write = 0
#else
      binary_write = 1
#endif

      call file_2b(ndn,'wallstress',3,1,tw_sfi,2,binary_write,0,1,ierr)
      call file_2b(ndn,'wallstress',3,1,tw_nfi,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',1,2,tw_lfi,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',1,2,tw_rfi,2,binary_write,1,1,ierr)

      call file_2b(ndn,'wallstress',3,2,tw_wfj,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',3,2,tw_efj,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',1,2,tw_lfj,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',1,2,tw_rfj,2,binary_write,1,1,ierr)

      call file_2b(ndn,'wallstress',3,2,tw_wfk,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',3,2,tw_efk,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',3,1,tw_sfk,2,binary_write,1,1,ierr)
      call file_2b(ndn,'wallstress',3,1,tw_nfk,2,binary_write,2,1,ierr)

      return
      end
