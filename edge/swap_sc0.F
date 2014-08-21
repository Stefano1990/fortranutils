C=======================================================================
C EXCHANGE INFORMATION FOR SCALAR VARIABLES ON PERIODIC BOUNDARIES AND
C BLOCK INTERFACES
C=======================================================================

      subroutine swap_sc(ivar1,igrid,var,work,mg_or_sg)

c     use mpi
      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "bc_type.inc"
#include "les.inc"
#include "data_type.inc"

C INCLUDE

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "comm.inc"
#include "cb.comm.inc"


C INPUT/OUTPUT

      integer ivar1, igrid
      real    var(*), work(lrwk2d)
      integer mg_or_sg

C LOCAL

      integer bc
      integer iblk, iside, iblk_next, iside_next
      integer d1, d2, d3, d1_next, d2_next, d3_next, p3
      integer side_mod, side_mod_next
      integer i1, i, j, k, io, jo, ko, ioff, ikoff
      integer ijko(3)
      integer nip1, njp1, nkp1, nip1_next, njp1_next, nkp1_next
      real    work2(lrwk2d)

      integer nsend, nrecv, ierr, iblk0
      integer ijkoff
      integer status(MPI_STATUS_SIZE)
      integer idisp2
#ifdef MPI1
      integer sreq(7*maxbl), i0
      integer idisp
#else
      integer(kind=mpi_address_kind) idisp
#endif

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      nsend = 0
      nrecv = 0
      iblk = 0

      do 10 iblk0 = 1,numbl_total

        if(pe.eq.bl2pe(iblk0)) then

          iblk = iblk + 1
          do 15 iside = 1,6

C-----------------------------------------------------------------------
C NOTE: ASSUME THAT GRID LINE ORIENTATION DIRECTION WILL NOT CHANGE
C   CROSS BLOCK BOUNDARIES, THAT IS, INCREASING DIRECTION WILL BE STILL
C   KEPT INCREASING. BUT THE ORIENTATION CAN CHANGE, FOR INSTANCE, I
C   GRID LINE CAN CHANGE TO J GRID LINE CROSS BLOCK BOUNDARIES
C-----------------------------------------------------------------------

            bc = bl_bc(iside,iblk)

            if (bc .eq. BCT_INTERFACE .or. bc .eq. BCT_PERIODIC) then

              iblk_next  = bl2bl(iside,iblk0)
              iside_next = bl2side(iside,iblk0)

              njp1 = iniknum(1,iside,iblk0,igrid)
              nkp1 = iniknum(2,iside,iblk0,igrid)

C-----------------------------------------------------------------------
C STORE THE DATA IN THE WORKSPACE
C-----------------------------------------------------------------------

              if (bc .eq. BCT_INTERFACE) then

                if (mg_or_sg .eq. 1) then
                  ioff=inikoff(1,iside,iblk0,igrid)
                else
                  ioff=inikoff(2,iside,iblk0,igrid)
                end if
                jo = inikoff(5,iside,iblk0,igrid)
                ko = inikoff(6,iside,iblk0,igrid)

                i1 = inwoff(iside,iblk0,igrid)
                do k = 1,nkp1
                  ikoff = ko*(k-1) + ioff
                  do j = 1,njp1
                    work(i1) = var(ikoff+jo*(j-1))
                    i1 = i1 + 1
                  end do
                end do

C PERIODIC B.C.: A CONSTANT (SUCH AS PRESSURE DROP) MAY BE INVOLVED

              else


                if (mg_or_sg .eq. 1) then
                  ioff=inikoff(1,iside,iblk0,igrid)
                else
                  ioff=inikoff(2,iside,iblk0,igrid)
                end if
                jo = inikoff(5,iside,iblk0,igrid)
                ko = inikoff(6,iside,iblk0,igrid)

                i1 = inwoff(iside,iblk0,igrid)

                do k = 1,nkp1
                  ikoff = ko*(k-1) + ioff
                  do j = 1,njp1
                    work(i1) = var(ikoff+jo*(j-1))
                    i1 = i1 + 1
                  end do
                end do

              end if
#ifdef MPI1
              if(bl2pe(iblk_next) .ne. pe) then
                nsend = nsend + 1
                i0 = inwoff(iside,iblk0,igrid)
                call MPI_ISEND(work(i0), i1-i0, REAL_TYPE_MPI,
     &          bl2pe(iblk_next), 7*(iblk_next-1)+iside_next,
     &          MPI_COMM_WORLD, sreq(nsend), ierr)
              endif
#endif

            end if

 15       continue

        end if

 10   continue

C-----------------------------------------------------------------------
C UPDATE INTERFACES OR PERIODIC BOUNDARIES
C-----------------------------------------------------------------------

#if defined SHMEM && ! defined DEBUG
      call shmem_barrier_all
#elif ! defined MPI1 && ! defined DEBUG
      call MPI_WIN_FENCE(0, workwin, ierr)
#endif


      iblk = 0
      do 20 iblk0 = 1,numbl_total

        if(pe.eq.bl2pe(iblk0)) iblk = iblk + 1

        do 25 iside = 1,6
#if defined SHMEM && defined DEBUG
          call shmem_barrier_all
#elif ! defined MPI1 && defined DEBUG
          call MPI_WIN_FENCE(0, workwin, ierr)
#endif

          iblk_next  = bl2bl(iside,iblk0)
          iside_next = bl2side(iside,iblk0)

          if(bl2bl(iside,iblk0).ne.0) then

            if(pe.eq.bl2pe(iblk0)) then

              bc = bl_bc(iside,iblk)

              if (bc .eq. BCT_INTERFACE .or. bc .eq. BCT_PERIODIC) then

                njp1 = iniknum(1,iside,iblk0,igrid)
                nkp1 = iniknum(2,iside,iblk0,igrid)

                if(pe.ne.bl2pe(iblk_next)) then
                  idisp = inwoff(iside_next,iblk_next,igrid)-1
                  idisp2 = inwoff(iside,iblk0,igrid)
#ifdef MPI1
                  call MPI_RECV(work2(idisp2), njp1*nkp1,
     &            REAL_TYPE_MPI, bl2pe(iblk_next), 7*(iblk0-1)+iside,
     &            MPI_COMM_WORLD, status, ierr)
#elif defined SHMEM
                  call shmem_get8(work2(idisp2), work(idisp+1),
     &            njp1*nkp1, bl2pe(iblk_next) )
#else
                  call MPI_GET(work2(idisp2), njp1*nkp1,
     &            REAL_TYPE_MPI, bl2pe(iblk_next), idisp, njp1*nkp1,
     &            REAL_TYPE_MPI, workwin, ierr)
#endif
                end if

              end if
            endif

          endif

#if defined SHMEM && defined DEBUG
          call shmem_barrier_all
#elif ! defined MPI1 && defined DEBUG
          call MPI_WIN_FENCE(0, workwin, ierr)
#endif
 25     continue
 20   continue



#if defined SHMEM && ! defined DEBUG
      call shmem_barrier_all
#elif ! defined MPI1 && ! defined DEBUG
      call MPI_WIN_FENCE(0, workwin, ierr)
#endif


      iblk = 0
      do 60 iblk0 = 1,numbl_total

        if(pe.eq.bl2pe(iblk0))iblk = iblk + 1

        do 50 iside = 1,6

          iblk_next  = bl2bl(iside,iblk0)
          iside_next = bl2side(iside,iblk0)

          if(pe.eq.bl2pe(iblk0))then

            bc = bl_bc(iside,iblk)
            if (bc .eq. BCT_INTERFACE .or. bc .eq. BCT_PERIODIC) then

              if(bl2bl(iside,iblk0).ne.0) then

                njp1 = iniknum(1,iside,iblk0,igrid)
                nkp1 = iniknum(2,iside,iblk0,igrid)
                if (mg_or_sg .eq. 1) then
                  ioff=inikoff(3,iside,iblk0,igrid)
                else
                  ioff=inikoff(4,iside,iblk0,igrid)
                end if
                jo = inikoff(5,iside,iblk0,igrid)
                ko = inikoff(6,iside,iblk0,igrid)

                if(pe.eq.bl2pe(iblk_next))then
                  i1 = inwoff(iside_next,iblk_next,igrid)
                  do k = 1,nkp1
                    ikoff = ko*(k-1) + ioff
                    do j = 1,njp1
                      var(ikoff+jo*(j-1)) = work(i1)
                      i1 = i1 + 1
                    end do
                  end do

                else
                  idisp2 = inwoff(iside,iblk0,igrid)
                  i1 = idisp2
                  do k = 1,nkp1
                    ikoff = ko*(k-1) + ioff
                    do j = 1,njp1
                      var(ikoff+jo*(j-1)) = work2(i1)
                      i1 = i1 + 1
                    end do
                  end do

                end if
              end if
            endif
          endif

 50     continue
 60   continue


#ifdef MPI1
      do i=1,nsend
        call MPI_WAIT(sreq(i), status, ierr)
      end do
#endif

      return
      end




