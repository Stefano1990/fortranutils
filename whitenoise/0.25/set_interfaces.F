C=======================================================================
C EXCHANGE INFORMATION FOR SCALAR VARIABLES ON PERIODIC BOUNDARIES AND
C BLOCK INTERFACES
C=======================================================================

      subroutine set_interfaces

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "bc_type.inc"
#include "les.inc"
#include "data_type.inc"

C INCLUDE

#include "cb.address.inc"
#include "cb.big.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "comm.inc"
#include "cb.comm.inc"


C LOCAL

      integer bc
      integer iblk, iside, iblk_next, iside_next
      integer d1, d2, d3, d1_next, d2_next, d3_next
      integer side_mod, side_mod_next
      integer io, jo, ko, p31, p32, i2
      integer ijko(3)
      integer i1, i, j, k, igrid, ibloff
      integer nip1, njp1, nkp1, nip1_next, njp1_next, nkp1_next
      integer itemp2(3),itemp(19*maxbl)

      integer nsend, nrecv, ierr, iblk0, iblk2
      integer sreq(6*maxbl), rreq(18*maxbl)
      integer istat(MPI_STATUS_SIZE)
#ifdef MPI1
      integer isize
#else
      integer(kind=MPI_ADDRESS_KIND) isize
#endif

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      call MPI_TYPE_EXTENT(REAL_TYPE_MPI, rlext, ierr)
      idwrcount = 0
      iccount = 0
      call MPI_ALLREDUCE(lrwk, ilrwk, 1, INTEGER_TYPE_MPI, MPI_MIN,
     &  MPI_COMM_WORLD, ierr)

      isize = rlext*ilrwk
#if ! defined MPI1 && ( ! defined SHMEM)
      call MPI_WIN_CREATE(ra(jw1a), isize, rlext,
     &  MPI_INFO_NULL, MPI_COMM_WORLD, workwin, ierr)
#endif

c     write(6,*)'ilrwk: ',ilrwk

      do iblk=1,numbl
        iblkl2g(iblk)=-1
      end do
      iblk = 0
      do iblk0 = 1,numbl_total
        if(pe.eq.bl2pe(iblk0))then
          iblk=iblk+1
          iblkg2l(iblk0)=iblk
          iblkl2g(iblk)=iblk0
        else
          iblkg2l(iblk0)=-1
        endif
      end do


      do igrid=1,nummgl
      i1 = 1
      iblk = 0
      nsend = 0
      nrecv = 0

      do 10 iblk0 = 1,numbl_total

        if(pe.eq.bl2pe(iblk0))iblk = iblk + 1
        do 15 iside = 1,6

C-----------------------------------------------------------------------
C NOTE: ASSUME THAT GRID LINE ORIENTATION DIRECTION WILL NOT CHANGE
C   CROSS BLOCK BOUNDARIES, THAT IS, INCREASING DIRECTION WILL BE STILL
C   KEPT INCREASING. BUT THE ORIENTATION CAN CHANGE, FOR INSTANCE, I
C   GRID LINE CAN CHANGE TO J GRID LINE CROSS BLOCK BOUNDARIES
C-----------------------------------------------------------------------

          iblk_next  = bl2bl(iside,iblk0)
          iside_next = bl2side(iside,iblk0)

          if(iblk_next.ne.0) then

          if(pe.eq.bl2pe(iblk0)) then

            bc = bl_bc(iside,iblk)

C-----------------------------------------------------------------------
C SINCE THE RESTRICTION ON THE CHANGE OF GRID LINE ORIENTATION, iside +
C isidec MUST BE ODD, I.E., iside AND isidec CANNOT BE BOTH ODD OR BOTH
C EVEN
C-----------------------------------------------------------------------

              side_mod      = mod(iside,2)
              side_mod_next = mod(iside_next,2)

              if (side_mod .eq. side_mod_next) then
                write(*,'(''*** Error: Interface coupling invalid: '',
     &                    ''bl '',i3, '' side '',i3,'' with bl '',i3,
     &                    ''side '',i3)')
     &                    iblk0, iside, iblk_next, iside_next
               stop
              end if

              d1 = (iside+1)/2
              d2 = cyclic(1,d1)
              d3 = cyclic(2,d1)

              nip1 = mg_nijkp1(d1,iblk,igrid)
              njp1 = mg_nijkp1(d2,iblk,igrid)
              nkp1 = mg_nijkp1(d3,iblk,igrid)

              ijko(1) = 1
              ijko(2) = mg_nijkp1(1,iblk,igrid)
              ijko(3) = mg_nijkp1(1,iblk,igrid)*
     &                  mg_nijkp1(2,iblk,igrid)

              io = ijko(d1)
              jo = ijko(d2)
              ko = ijko(d3)

              p31 = mg_off(iblk,igrid)
              p32 = sg_off(iblk)

              if (side_mod .eq. 0) then
                i = nip1 - 1
                i2 = nip1
              else
                i = 2
                i2 = 1
              end if

              inikoff(1,iside,iblk0,igrid) = 1 + io*(i-1) + p31
              inikoff(2,iside,iblk0,igrid) = 1 + io*(i-1) + p32
              inikoff(3,iside,iblk0,igrid) = 1 + io*(i2-1) + p31
              inikoff(4,iside,iblk0,igrid) = 1 + io*(i2-1) + p32
              inikoff(5,iside,iblk0,igrid) = jo
              inikoff(6,iside,iblk0,igrid) = ko

              if(pe.eq.bl2pe(iblk_next))then
                d1_next = (iside_next+1)/2
                d2_next = cyclic(1,d1_next)
                d3_next = cyclic(2,d1_next)
                nip1_next = mg_nijkp1(d1_next,iblk,igrid)
                njp1_next = mg_nijkp1(d2_next,iblk,igrid)
                nkp1_next = mg_nijkp1(d3_next,iblk,igrid)
              else
                CALL MPI_RECV(itemp2, 3, INTEGER_TYPE_MPI,
     &            bl2pe(iblk_next), iblk0*6 + iside, MPI_COMM_WORLD,
     &            istat, ierr)
                nip1_next = itemp2(1)
                njp1_next = itemp2(2)
                nkp1_next = itemp2(3)
              endif

C CHECK IF THE TWO BLOCKS MATCH
C
              if (njp1.ne.njp1_next .or. nkp1.ne.nkp1_next) then
                write(*,'(''*** Error: Block '',i3,'' face '',i3,
     &            '' to Block '',i3,'' face '',i3,'' mismatch'')')
     &            iblk0,iside,iblk_next,iside_next
              end if

              inwoff(iside,iblk0,igrid) = i1
              iniknum(1,iside,iblk0,igrid) = njp1
              iniknum(2,iside,iblk0,igrid) = nkp1


              if(pe.ne.bl2pe(iblk_next))then
                call MPI_SEND(inwoff(iside,iblk0,igrid),
     &          1, INTEGER_TYPE_MPI, bl2pe(iblk_next), iblk0*6+iside,
     &          MPI_COMM_WORLD, ierr)
                call MPI_SEND(inikoff(1,iside,iblk0,igrid),
     &          3, INTEGER_TYPE_MPI, bl2pe(iblk_next), iblk0*6+iside,
     &          MPI_COMM_WORLD, ierr)
                call MPI_SEND(iniknum(1,iside,iblk0,igrid),
     &          2, INTEGER_TYPE_MPI, bl2pe(iblk_next), iblk0*6+iside,
     &          MPI_COMM_WORLD, ierr)
              endif

              i1=i1+njp1*nkp1+1

          elseif(pe.eq.bl2pe(iblk_next)) then

            iblk2=iblkg2l(iblk_next)
            bc = bl_bc(iside_next,iblk2)

              d1_next = (iside_next+1)/2
              d2_next = cyclic(1,d1_next)
              d3_next = cyclic(2,d1_next)
              ibloff = (iblk2-1)*18+3*(iside-1)+1
              itemp(ibloff) = mg_nijkp1(d1_next,iblk2,igrid)
              itemp(ibloff+1) = mg_nijkp1(d2_next,iblk2,igrid)
              itemp(ibloff+2) = mg_nijkp1(d3_next,iblk2,igrid)

              nsend = nsend + 1
              call MPI_ISEND(itemp(ibloff), 3, INTEGER_TYPE_MPI,
     &        bl2pe(iblk0), iblk0*6+iside, MPI_COMM_WORLD,
     &        sreq(nsend), ierr)

              nrecv = nrecv + 1
              call MPI_IRECV(inwoff(iside,iblk0,igrid),
     &        1, INTEGER_TYPE_MPI, bl2pe(iblk0), iblk0*6+iside,
     &        MPI_COMM_WORLD, rreq(nrecv), ierr)
              nrecv = nrecv + 1
              call MPI_IRECV(inikoff(1,iside,iblk0,igrid),
     &        3, INTEGER_TYPE_MPI, bl2pe(iblk0), iblk0*6+iside,
     &        MPI_COMM_WORLD, rreq(nrecv), ierr)
              nrecv = nrecv + 1
              call MPI_IRECV(iniknum(1,iside,iblk0,igrid),
     &        2, INTEGER_TYPE_MPI, bl2pe(iblk0), iblk0*6+iside,
     &        MPI_COMM_WORLD, rreq(nrecv), ierr)

          end if
          end if

 15     continue

 10   continue

      inoffmax(igrid) = i1+1

      do i=1,nsend
        call MPI_WAIT(sreq(i), istat, ierr)
      end do

      do i=1,nrecv
        call MPI_WAIT(rreq(i), istat, ierr)
      end do

      end do

      return
      end





