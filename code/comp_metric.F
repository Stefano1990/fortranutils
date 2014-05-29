C=======================================================================
C METRIC QUANTITY CALCULATION CONTROL
C=======================================================================

      subroutine comp_metric(x,y,t,
     &  dip,diu,div,dix,rdip,rdiv,rdix,
     &  djp,dju,djv,djx,rdjp,rdju,rdjx,
     &  hijp,hjip,hkip,hkjp,
     &  hiju,hjiu,hkiu,hkju,
     &  hijv,hjiv,hkiv,hkjv,
     &  rp,ru,rv,
     &  didjp,rdidjp,didju,rdidju,djdiv,rdjdiv,
     &  dkp,dkw,
     &  swap_geo,pnti)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

C INCLUDE
#include "cb.gridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

C INPUT/OUTPUT PARAMETERS

      real    x(*), y(*), t(*)
      real    dip(*), diu(*), div(*), dix(*), rdip(*), rdiv(*), rdix(*)
      real    djp(*), dju(*), djv(*), djx(*), rdjp(*), rdju(*), rdjx(*)
      real    hijp(*), hjip(*), hkip(*), hkjp(*)
      real    hiju(*), hjiu(*), hkiu(*), hkju(*)
      real    hijv(*), hjiv(*), hkiv(*), hkjv(*)
      real    rp(*), ru(*), rv(*)
      real    didjp(*), rdidjp(*), didju(*), rdidju(*),
     &        djdiv(*), rdjdiv(*)
      real    dkp(*), dkw(*)
      real    swap_geo(*)


C LOCAL QUANTITIES

      integer iblk, ilev, iside, iblk_next, iside_next, ierr
      integer nip1, njp1, nkp1, nip2, njp2, nkp2
      integer p_ij, p_k, pnt
      integer p_xy, p_z

      integer status(MPI_STATUS_SIZE), sreq(18*(maxbl+1)), iblk0
c     integer tag(18*(maxbl+1))

      integer k, nsend
      real    swap_geo2(34*maxijk1)
c     logical lprobe


C INPUT/OUTPUT PARAMETERS
      integer pnti(6,numbl_total)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do ilev = 1,nummgl

C-----------------------------------------------------------------------
C COMPUTE THE METRIC QUANTITIES FOR ALL BLOCKS ON THE SAME GRID LEVEL
C-----------------------------------------------------------------------

        nsend = 0
        pnt = 1
        iblk = 0

        do 15 iblk0 = 1,numbl_total

          if(pe.eq.bl2pe(iblk0))then
            iblk = iblk + 1

            p_ij = mg_off_2d_ij(iblk,ilev) + 1
            p_k  = mg_off_1d_k(iblk,ilev)  + 1
            p_xy = mg_off_2d_xy(iblk,ilev) + 1
            p_z  = mg_off_1d_z(iblk,ilev)  + 1

            nip1 = mg_nijkp1(1,iblk,ilev)
            njp1 = mg_nijkp1(2,iblk,ilev)
            nkp1 = mg_nijkp1(3,iblk,ilev)
            nip2 = nip1 + 1
            njp2 = njp1 + 1
            nkp2 = nkp1 + 1

C-----------------------------------------------------------------------
C CALCULATE METRIC QUANTITIES
C-----------------------------------------------------------------------

            call metric_block(nip1,njp1,nkp1,nip2,njp2,nkp2,
     &         x(p_xy),y(p_xy),t(p_z),
     &         dip(p_ij),diu(p_ij),div(p_ij),dix(p_ij),rdip(p_ij),
     &         rdiv(p_ij),rdix(p_ij),
     &         djp(p_ij),dju(p_ij),djv(p_ij),djx(p_ij),rdjp(p_ij),
     &         rdju(p_ij),rdjx(p_ij),
     &         hijp(p_ij),hjip(p_ij),hkip(p_ij),hkjp(p_ij),
     &         hiju(p_ij),hjiu(p_ij),hkiu(p_ij),hkju(p_ij),
     &         hijv(p_ij),hjiv(p_ij),hkiv(p_ij),hkjv(p_ij),
     &         rp(p_ij),ru(p_ij),rv(p_ij),
     &         didjp(p_ij),rdidjp(p_ij),didju(p_ij),rdidju(p_ij),
     &         djdiv(p_ij),rdjdiv(p_ij),
     &         dkp(p_k),dkw(p_k),gotrad,iblk0)

C-----------------------------------------------------------------------
C STORE THE VALUES ON PERIODIC BOUNDARIES OR INTERFACES IN WORKSPACE
C-----------------------------------------------------------------------

            do 20 iside = 1,6

c             goto 20
  
              iblk_next  = bl2bl(iside,iblk0)
              iside_next = bl2side(iside,iblk0)

              if (iblk_next.gt.0) then
                pnti(iside_next,iblk_next) = pnt 

C SWAP FOR METRIC QUANTITIES NEEDS SPECIAL TREATMENT BECAUSE THERE ARE
C ONE- AND TWO-DIMENSIONAL ARRAYS

                call metric_swap1(iside,bl_bc(iside,iblk),nip1,njp1,
     &             nkp1,dip(p_ij),diu(p_ij),div(p_ij),dix(p_ij),
     &             rdip(p_ij),rdiv(p_ij),rdix(p_ij),
     &             djp(p_ij),dju(p_ij),djv(p_ij),djx(p_ij),
     &             rdjp(p_ij),rdju(p_ij),rdjx(p_ij),
     &             hijp(p_ij),hjip(p_ij),hkip(p_ij),hkjp(p_ij),
     &             hiju(p_ij),hjiu(p_ij),hkiu(p_ij),hkju(p_ij),
     &             hijv(p_ij),hjiv(p_ij),hkiv(p_ij),hkjv(p_ij),
     &             rp(p_ij),ru(p_ij),rv(p_ij),
     &             didjp(p_ij),didju(p_ij),rdidju(p_ij),djdiv(p_ij),
     &             rdjdiv(p_ij),
     &             dkp(p_k),dkw(p_k),
     &             pnt,swap_geo,iblk_next)

                if(bl2pe(iblk_next).ne.pe)then

                  nsend = nsend + 1
                  call MPI_ISEND(pnt-pnti(iside_next,iblk_next), 1,
     &            INTEGER_TYPE_MPI, bl2pe(iblk_next),
     &            13*maxbl*(ilev-1)+13*(iblk0-1)+2*iside-1,
     &            MPI_COMM_WORLD, sreq(nsend), ierr)

                  nsend = nsend + 1
                  call MPI_ISEND(swap_geo(pnti(iside_next,iblk_next)),
     &            pnt-pnti(iside_next,iblk_next), REAL_TYPE_MPI,
     &            bl2pe(iblk_next),
     &            13*maxbl*(ilev-1)+13*(iblk0-1)+2*iside,
     &            MPI_COMM_WORLD, sreq(nsend), ierr)
                end if

              end if

 20         continue

          end if

 15     continue



C-----------------------------------------------------------------------
C EXCHANGE THE VALUES FOR PERIODIC OR INTERFACES
C-----------------------------------------------------------------------

        iblk = 0
        do 25 iblk0 = 1,numbl_total
  
          if(bl2pe(iblk0).eq.pe) then
            iblk = iblk + 1

            p_ij = mg_off_2d_ij(iblk,ilev) + 1
            p_k  = mg_off_1d_k(iblk,ilev)  + 1

            nip1 = mg_nijkp1(1,iblk,ilev)
            njp1 = mg_nijkp1(2,iblk,ilev)
            nkp1 = mg_nijkp1(3,iblk,ilev)

            do 30 iside = 1,6

              iblk_next  = bl2bl(iside,iblk0)
              iside_next = bl2side(iside,iblk0)

              if (iblk_next.gt.0) then

                if(bl2pe(iblk_next).ne.pe)then

                  call MPI_RECV(pnt, 1, INTEGER_TYPE_MPI,
     &              bl2pe(iblk_next),
     &              13*maxbl*(ilev-1)+13*(iblk_next-1)+2*iside_next-1,
     &              MPI_COMM_WORLD, status, ierr)

                  call MPI_RECV(swap_geo2, pnt, REAL_TYPE_MPI,
     &              bl2pe(iblk_next),
     &              13*maxbl*(ilev-1)+13*(iblk_next-1)+2*iside_next,
     &              MPI_COMM_WORLD, status, ierr)
                  pnt = 1
                  call metric_swap2(iside,bl_bc(iside,iblk),nip1,njp1,
     &              nkp1,dip(p_ij),diu(p_ij),div(p_ij),dix(p_ij),
     &              rdip(p_ij),rdiv(p_ij),rdix(p_ij),
     &              djp(p_ij),dju(p_ij),djv(p_ij),djx(p_ij),
     &              rdjp(p_ij),rdju(p_ij),rdjx(p_ij),
     &              hijp(p_ij),hjip(p_ij),hkip(p_ij),hkjp(p_ij),
     &              hiju(p_ij),hjiu(p_ij),hkiu(p_ij),hkju(p_ij),
     &              hijv(p_ij),hjiv(p_ij),hkiv(p_ij),hkjv(p_ij),
     &              rp(p_ij),ru(p_ij),rv(p_ij),
     &              didjp(p_ij),didju(p_ij),rdidju(p_ij),djdiv(p_ij),
     &              rdjdiv(p_ij),
     &              dkp(p_k),dkw(p_k),
     &              pnt,swap_geo2,iblk0)

                else
                  pnt = pnti(iside,iblk0)
                  call metric_swap2(iside,bl_bc(iside,iblk),nip1,njp1,
     &              nkp1,dip(p_ij),diu(p_ij),div(p_ij),dix(p_ij),
     &              rdip(p_ij),rdiv(p_ij),rdix(p_ij),
     &              djp(p_ij),dju(p_ij),djv(p_ij),djx(p_ij),
     &              rdjp(p_ij),rdju(p_ij),rdjx(p_ij),
     &              hijp(p_ij),hjip(p_ij),hkip(p_ij),hkjp(p_ij),
     &              hiju(p_ij),hjiu(p_ij),hkiu(p_ij),hkju(p_ij),
     &              hijv(p_ij),hjiv(p_ij),hkiv(p_ij),hkjv(p_ij),
     &              rp(p_ij),ru(p_ij),rv(p_ij),
     &              didjp(p_ij),didju(p_ij),rdidju(p_ij),djdiv(p_ij),
     &              rdjdiv(p_ij),
     &              dkp(p_k),dkw(p_k),
     &              pnt,swap_geo,iblk0)

                end if

              end if

 30         continue

          end if

 25     continue



        do k=1,nsend
          call MPI_WAIT(sreq(k), status, ierr)
        end do

      enddo  !ilev

      return
      end
