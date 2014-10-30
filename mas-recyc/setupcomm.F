C=======================================================================
C COMMUNICATIONS INITIALIZATION ROUTINE
C=======================================================================

      subroutine setupcomm()
 
C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "comm.inc"

C LOCAL
        integer ierr

C COMMON
 
#include "cb.comm.inc"




      call MPI_INIT(ierr)

      call MPI_COMM_RANK(MPI_COMM_WORLD, pe, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numpe, ierr)

      if(pe.eq.0)master = .true.





      return
      end

