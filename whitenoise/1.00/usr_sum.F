
C Routines:

c usr_sum_double  - Use REAL*8 throughout
c usr_sum_single  - Use REAL*4 throughout
c usr_sum_single2 - Input and output REAL*4, intermediate arithmetic REAL*8
c usr_sum_double2 - Input and output REAL*8, intermediate arithmetic REAL*16

 
c All the routines have the same arguments, which are

c sendbuf      - array of size numbl containing terms
c recvbuf      - scalar floating point to hold result
c count        - scalar integer equal to numbl_total
c datatype     - MPI datatype in sendbuf and recvbuf
c communicator - MPI communicator
c root         - rank of root process for this operation
c order        - bl2pe(numbl_total)
c rank         - rank of calling process
c broadcast    - .true. => make this an MPI_ALLREDUCE, else not
c ierr         - in case I need to return status info



       subroutine usr_sum_double(sendbuf, recvbuf, count, datatype,
     &                           communicator, root, order, rank,
     &                           broadcast, ierr)

       implicit none

#include "mpif.h"

       logical broadcast
       integer count, datatype, communicator, root, rank, ierr, i
       integer order(count), request(count)
       integer icount

       integer, parameter :: kind1 = selected_real_kind(13,300)
       integer, parameter :: kind2 = selected_real_kind(13,300)

       real(kind=kind1) :: sendbuf(*), recvbuf
       real(kind=kind2) :: interm

       integer status(MPI_STATUS_SIZE)

       icount = 0
#ifdef T3E
       interm = 0.0
#else
       interm = 0.0D0
#endif

       do i=1,count
         if (order(i) .eq. rank) then
            icount = icount + 1
            call MPI_ISEND(sendbuf(icount), 1, datatype, root, i,
     &                     communicator, request(icount), ierr)
         end if
         if (root .eq. rank) then
            call MPI_RECV(recvbuf, 1, datatype, order(i), i,
     &                    communicator, status, ierr)
            interm = interm + recvbuf
         end if
       end do

       do i=1,icount
         call MPI_WAIT(request(i), status, ierr)
       end do

       if(rank .eq. root) recvbuf = interm

       if(broadcast) then
         call MPI_BCAST(recvbuf, 1, datatype, root, communicator, ierr)
       end if

       return
       end







       subroutine usr_sum_single(sendbuf, recvbuf, count, datatype,
     &                           communicator, root, order, rank,
     &                           broadcast, ierr)

       implicit none

#include "mpif.h"

       logical broadcast
       integer count, datatype, communicator, root, rank, ierr, i
       integer order(count), request(count)
       integer icount

       integer, parameter :: kind1 = selected_real_kind(6,30)
       integer, parameter :: kind2 = selected_real_kind(6,30)

       real(kind=kind1) :: sendbuf(*), recvbuf
       real(kind=kind2) :: interm

       integer status(MPI_STATUS_SIZE)

       icount = 0
#ifdef T3E
       interm = 0.0
#else
       interm = 0.0D0
#endif

       do i=1,count
         if (order(i) .eq. rank) then
            icount = icount + 1
            call MPI_ISEND(sendbuf(icount), 1, datatype, root, i,
     &                     communicator, request(icount), ierr)
         end if
         if (root .eq. rank) then
            call MPI_RECV(recvbuf, 1, datatype, order(i), i,
     &                    communicator, status, ierr)
            interm = interm + recvbuf
         end if
       end do

       do i=1,icount
         call MPI_WAIT(request(i), status, ierr)
       end do

       if(rank .eq. root) recvbuf = interm

       if(broadcast) then
         call MPI_BCAST(recvbuf, 1, datatype, root, communicator, ierr)
       end if

       return
       end







       subroutine usr_sum_single2(sendbuf, recvbuf, count, datatype,
     &                            communicator, root, order, rank,
     &                            broadcast, ierr)

       implicit none

#include "mpif.h"

       logical broadcast
       integer count, datatype, communicator, root, rank, ierr, i
       integer order(count), request(count)
       integer icount

       integer, parameter :: kind1 = selected_real_kind(6,30)
       integer, parameter :: kind2 = selected_real_kind(13,300)

       real(kind=kind1) :: sendbuf(*), recvbuf
       real(kind=kind2) :: interm

       integer status(MPI_STATUS_SIZE)

       icount = 0
#ifdef T3E
       interm = 0.0
#else
       interm = 0.0D0
#endif

       do i=1,count
         if (order(i) .eq. rank) then
            icount = icount + 1
            call MPI_ISEND(sendbuf(icount), 1, datatype, root, i,
     &                     communicator, request(icount), ierr)
         end if
         if (root .eq. rank) then
            call MPI_RECV(recvbuf, 1, datatype, order(i), i,
     &                    communicator, status, ierr)
            interm = interm + recvbuf
         end if
       end do

       do i=1,icount
         call MPI_WAIT(request(i), status, ierr)
       end do

       if(rank .eq. root) recvbuf = interm

       if(broadcast) then
         call MPI_BCAST(recvbuf, 1, datatype, root, communicator, ierr)
       end if

       return
       end







       subroutine usr_sum_double2(sendbuf, recvbuf, count, datatype,
     &                            communicator, root, order, rank,
     &                            broadcast, ierr)

       implicit none

#include "mpif.h"

       logical broadcast
       integer count, datatype, communicator, root, rank, ierr, i
       integer order(count), request(count)
       integer icount

       integer, parameter :: kind1 = selected_real_kind(13,300)
       integer, parameter :: kind2 = selected_real_kind(20,400)

       real(kind=kind1) :: sendbuf(*), recvbuf
#if defined DEC
       real(kind=kind2) :: interm
#elif defined T3E
       real(kind=8)     :: interm
#elif defined SGI
       real(kind=16)    :: interm
#endif

       integer status(MPI_STATUS_SIZE)

       icount = 0
#if defined T3E
       interm = 0.0
#else
       interm = 0.0D0
#endif

       do i=1,count
         if (order(i) .eq. rank) then
            icount = icount + 1
            call MPI_ISEND(sendbuf(icount), 1, datatype, root, i,
     &                     communicator, request(icount), ierr)
         end if
         if (root .eq. rank) then
            call MPI_RECV(recvbuf, 1, datatype, order(i), i,
     &                    communicator, status, ierr)
            interm = interm + recvbuf
         end if
       end do

       do i=1,icount
         call MPI_WAIT(request(i), status, ierr)
       end do

       if(rank .eq. root) recvbuf = interm

       if(broadcast) then
         call MPI_BCAST(recvbuf, 1, datatype, root, communicator, ierr)
       end if

       return
       end






