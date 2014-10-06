C=======================================================================
C EXCHANGE x, y, z COORDINATES ON PERIODIC BOUNDARIES AND BLOCK
C INTERFACES
C=======================================================================

C NOTE: THIS ROUTINE IS NOT GENERAL. THE ORIENTATION OF GRID LINES ARE
C   NOT ALLOWED TO CHANGE CROSS BLOCK BOUNDARIES. IT IS NOT IMPLEMENTED
C   IN A GENERAL WAY BECAUSE THE PACE PART IS NOT GENERAL.

C-----------------------------------------------------------------------

      subroutine swap_xyz(coorname,ilev,x,work)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "bc_type.inc"
#include "data_type.inc"


C INPUT/OUTPUT

      integer   ilev
      real      x(*), work(*)
      character coorname

C LOCAL

      integer bc
      integer iblk, iside, iblk_next, iside_next
      integer d1, d1_next, d2, d2_next, pnt
      integer side_mod, side_mod_next
      integer nip1, njp1, nip2, njp2, nip1_next, njp1_next
      integer i1, i, j, ioff, joff
      integer maxijk1_coor

      integer ierr, nrecs, nsends, i0, iblk_local
      integer rreq(6*maxbl), sreq(6*maxbl), tag(6*maxbl)
      integer status(MPI_STATUS_SIZE)
      logical pebound

C INCLUDE

#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      i1 = 0
      nsends = 0
      iblk_local = 0
      maxijk1_coor = maxijk1 + 1

C-----------------------------------------------------------------------
C PUSH THE INFORMATION FOR EXCHANGE INTO WORKSPACE (ACTUALLY NOT NECE-
C SSARY ON SHARED MEMORY MACHINES)
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl_total

        if (pe.eq.bl2pe(iblk)) then

        iblk_local = iblk_local + 1

        do 15 iside = 1,6

C NOTE: ASSUME THAT GRID LINE ORIENTATION DIRECTION WILL NOT CHANGE
C   CROSS BLOCK BOUNDARIES, THAT IS, INCREASING DIRECTION WILL BE STILL
C   KEPT INCREASING. BUT THE ORIENTATION CAN CHANGE, FOR INSTANCE, I
C   GRID LINE CAN CHANGE TO J GRID LINE CROSS BLOCK BOUNDARIES

          bc = bl_bc(iside,iblk_local)

          if (bc .eq. BCT_INTERFACE) then

            iblk_next  = bl2bl(iside,iblk)
            iside_next = bl2side(iside,iblk)

            if(bl2pe(iblk).ne.bl2pe(iblk_next))then
              pebound = .true.
            else
              pebound = .false.
            endif

            if (iblk_next .gt. 0) then

C SINCE THE RESTRICTION ON THE CHANGE OF GRID LINE ORIENTATION, iside +
C iside_next MUST BE ODD, I.E., iside AND iside_next CANNOT BE BOTH ODD
C OR BOTH EVEN

              side_mod      = mod(iside     ,2)
              side_mod_next = mod(iside_next,2)

              if (side_mod .eq. side_mod_next) then
                write(*,'(''*** Error: Interface coupling invalid: '',
     &                    ''bl '',i3, '' side '',i3,'' with bl '',i3,
     &                    ''side '',i3)')
     &                    iblk, iside, iblk_next, iside_next
                stop
              end if

C THIS IS AD HOC (WHEN d1 = 3, d2 CAN BE SET TO ANY VALUE IN 1,2,3)

              d1 = (iside+1)/2
              if (d1 .eq. 1) then
                d2 = 2
              else if (d1.eq.2) then
                d2 = 1
              else
                d2 = 1
              end if
              nip1 = mg_nijkp1(d1,iblk_local,ilev)
              njp1 = mg_nijkp1(d2,iblk_local,ilev)

              d1_next = d1
              d2_next = d2
c             nip1_next = mg_nijkp1(d1_next,iblk_next,ilev)
c             njp1_next = mg_nijkp1(d2_next,iblk_next,ilev)

C CHECK IF THE TWO BLOCKS MATCH (ASSUME THAT IN DIRECTION 3 THE VOLUME
C IS FORMED BY EITHER PROJECTION OR RATATION)
C
C RJ: Left this out. Might put back later. Would require
C     communication of njp1_next.
C
c             if (njp1.ne.njp1_next) then
c               write(*,'(''*** Error: Block '',i3,'' face '',i3,
c    &            '' to Block '',i3,'' face '',i3,'' mismatch'')')
c    &            iblk,iside,iblk_next,iside_next
c               stop
c             end if

C - FOR THE x, y COORDINATES

              nip2 = nip1 + 1
              njp2 = njp1 + 1

              if ((coorname.eq.'x' .or. coorname.eq.'y') .and.
     &             iside.le.4) then

                pnt = mg_off_2d_xy(iblk_local,ilev)

c               i1 = 4*maxijk1_coor*(iblk_local-1)
                i1 = 4*maxijk1_coor*(iblk-1)
     &             + maxijk1_coor*(iside-1)
                i0 = i1

                if (d1 .eq. 2) then
                  if (side_mod .eq. 0) then
                    joff = pnt + njp2*(nip2-3)
                  else
                    joff = pnt + njp2*(3-1)
                  end if
                  do i = 1,njp2
                    i1 = i1 + 1
                    work(i1) = x(i+joff)
                  end do
                else
                  if (side_mod .eq. 0) then
                    ioff = pnt + nip2-2
                  else
                    ioff = pnt + 3
                  end if
                  do j = 1,njp2
                    i1 = i1 + 1
                    work(i1) = x(ioff+nip2*(j-1))
                  end do
                end if

C - FOR THE z COORDINATE

              else if(iside.gt.4) then

                pnt = mg_off_1d_z(iblk_local,ilev)

                i1 = 4*(iblk-1) + side_mod + 1
                i0 = i1 - 1

                if (side_mod .eq. 0) then
                  work(i1) = x(pnt+nip2-2)
                else
                  work(i1) = x(pnt+1)
                end if

              else
                pebound = .false.
              end if

C - CHECK USE OF WORKSPACE

              if (i1 .gt. lrwk) then
                write(*,'(''*** Error: '',
     &            ''workspace operation exceeding its length'')')
                stop
              end if

              if (pebound) then
                nsends = nsends + 1
                sreq(nsends) = 1024*iblk + iside
                tag(nsends) = 1024*iblk_next + iside_next
                call MPI_ISEND(work(i0+1),i1 - i0, REAL_TYPE_MPI,
     &            bl2pe(iblk_next), tag(nsends), MPI_COMM_WORLD,
     &            sreq(nsends), ierr)
              endif

            end if
          end if

 15     continue

        end if

 10   continue

C-----------------------------------------------------------------------
C POP THE INFORMATION FOR EXCHANGE OUT OF WORKSPACE AND WRITE INTO
C THE CORRESPONDING COORDINATE ARRAY
C-----------------------------------------------------------------------

      iblk_local = 0
      nrecs = 0

      do 20 iblk = 1,numbl_total

        if (pe.eq.bl2pe(iblk)) then

        iblk_local = iblk_local + 1

        do 25 iside = 1,6

          bc = bl_bc(iside,iblk_local)

          if (bc .eq. BCT_INTERFACE) then

            iblk_next  = bl2bl(iside,iblk)
            iside_next = bl2side(iside,iblk)

            if(bl2pe(iblk).ne.bl2pe(iblk_next))then
              pebound = .true.
            else
              pebound = .false.
            endif

            if (iblk_next .gt. 0) then

              side_mod      = mod(iside     ,2)
              side_mod_next = mod(iside_next,2)

              d1 = (iside+1)/2
              if (d1 .eq. 1) then
                d2 = 2
              else if (d1.eq.2) then
                d2 = 1
              else
                d2 = 1
              end if
              nip1 = mg_nijkp1(d1,iblk_local,ilev)
              njp1 = mg_nijkp1(d2,iblk_local,ilev)

              d1_next = d1
              d2_next = d2

              nip2 = nip1 + 1
              njp2 = njp1 + 1

C - FOR THE x, y COORDINATES

              if ((coorname.eq.'x' .or. coorname.eq.'y') .and.
     &             iside.le.4) then

                pnt = mg_off_2d_xy(iblk_local,ilev)

                if(pebound)then
c                 i1 = 0
                  i1 = 4*maxijk1_coor*(iblk_next-1) +
     &                   maxijk1_coor*(iside_next-1)
                  nrecs = nrecs + 1
                  tag(nrecs) = 1024*iblk + iside
                  rreq(nrecs) = 1024*iblk + iside
                  call MPI_RECV(work(i1+1), njp2, REAL_TYPE_MPI,
     &              bl2pe(iblk_next), tag(nrecs), MPI_COMM_WORLD,
     &              rreq(nrecs), ierr)
                else
                  i1 = 4*maxijk1_coor*(iblk_next-1) +
     &                   maxijk1_coor*(iside_next-1)
                endif

                if (d1 .eq. 2) then
                  if (side_mod .eq. 0) then
                    joff = pnt + njp2*(nip2-1)
                  else
                    joff = pnt
                  end if
                  do i = 1,njp2
                    i1 = i1 + 1
                    x(i+joff) = work(i1)
                  end do
                else
                  if (side_mod .eq. 0) then
                    ioff = pnt + nip2
                  else
                    ioff = pnt + 1
                  end if
                  do j = 1,njp2
                    i1 = i1 + 1
                    x(ioff+nip2*(j-1)) = work(i1)
                  end do
                end if

C - FOR THE z COORDINATE

              else if (iside.gt.4) then

                pnt = mg_off_1d_z(iblk_local,ilev)

                if(pebound)then
c                 i1 = 1
                  i1 = 4*(iblk_next-1) + side_mod_next + 1
                  nrecs = nrecs + 1
                  tag(nrecs) = 1024*iblk + iside
                  rreq(nrecs) = 1024*iblk + iside
                  call MPI_RECV(work(i1), 1, REAL_TYPE_MPI,
     &              bl2pe(iblk_next), tag(nrecs), MPI_COMM_WORLD,
     &              status, ierr)
                else
                  i1 = 4*(iblk_next-1) + side_mod_next + 1
                endif

                if (side_mod .eq. 0) then
                  x(pnt+nip2) = work(i1)
                else
                  x(pnt+1   ) = work(i1)
                end if

              else
                pebound = .false.
              end if

            end if
          end if

 25     continue

        end if

 20   continue

        do i = 1, nsends
          call MPI_WAIT(sreq(i), status, ierr)
        end do



      return
      end
