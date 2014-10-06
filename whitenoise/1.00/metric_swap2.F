C=======================================================================
C EXCHANGE METRIC INFORMATION FOR INTERFACES AND PERIODIC BOUNDARIES
C - 1. THE CURRENT CODING REQUIRES THAT THE GRID ORIENTATION IS THE SAME
C      CROSS BLOCK BOUNDARIES
C - 2. COPY THE INFORMATION IN THE WORKSPACE TO THE PROPER DESTINATION
C - 3. SWAPPING FOR CERTAIN QUANTITIES ARE NOT NECESSARY; SIMPLY SWAP-
C      PING ALL
C=======================================================================

      subroutine metric_swap2(iside,bc,nip1,njp1,nkp1,
     &  dip,diu,div,dix,rdip,rdiv,rdix,
     &  djp,dju,djv,djx,rdjp,rdju,rdjx,
     &  hijp,hjip,hkip,hkjp,
     &  hiju,hjiu,hkiu,hkju,
     &  hijv,hjiv,hkiv,hkjv,
     &  rp,ru,rv,
     &  didjp,didju,rdidju,djdiv,rdjdiv,
     &  dkp,dkw,
     &  pnti,swap_geo,iblk)

      implicit none

C PARAMETERS

#include "dim.inc"
#include "bc_type.inc"

C INPUT/OUTPUT

      integer  iside, pnti, bc, iblk
      integer  nip1, njp1, nkp1
      real     dip(nip1*njp1), diu(nip1*njp1), div(nip1*njp1),
     &         dix(nip1*njp1), rdip(nip1*njp1), rdiv(nip1*njp1),
     &         rdix(nip1*njp1)
      real     djp(nip1*njp1), dju(nip1*njp1), djv(nip1*njp1),
     &         djx(nip1*njp1), rdjp(nip1*njp1), rdju(nip1*njp1),
     &         rdjx(nip1*njp1)
      real     hijp(nip1*njp1), hjip(nip1*njp1), hkip(nip1*njp1),
     &         hkjp(nip1*njp1)
      real     hiju(nip1*njp1), hjiu(nip1*njp1), hkiu(nip1*njp1),
     &         hkju(nip1*njp1)
      real     hijv(nip1*njp1), hjiv(nip1*njp1), hkiv(nip1*njp1),
     &         hkjv(nip1*njp1)
      real     rp(nip1*njp1), ru(nip1*njp1), rv(nip1*njp1)
      real     didjp(nip1*njp1), didju(nip1*njp1), rdidju(nip1*njp1)
      real     djdiv(nip1*njp1), rdjdiv(nip1*njp1)
      real     dkp(nkp1), dkw(nkp1)
      real     swap_geo(*)

C LOCAL

      integer  i, j, ij, ioff, io, jo, je

C INCLUDE

#include "cb.lgrid.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      
C-----------------------------------------------------------------------
C MOVE INFORMATION OUT OF THE WORKSPACE
C-----------------------------------------------------------------------

      if (bc.eq.BCT_INTERFACE .or. bc.eq.BCT_PERIODIC) THEN

       if (iside .eq. 5) then

          dkp(1   ) = swap_geo(pnti  )
          dkw(1   ) = swap_geo(pnti+1)
          pnti = pnti + 2

       else if (iside .eq. 6) then

          dkp(nkp1)  = swap_geo(pnti  )
          dkw(nkp1)  = swap_geo(pnti+1)
          pnti = pnti + 2

        else

          if (iside .eq. 1) then
            i  = 1
            io = 1
            jo = nip1
            je = njp1
          else if (iside .eq. 2) then
            i  = nip1
            io = 1
            jo = nip1
            je = njp1
          else if (iside .eq. 3) then
            i  = 1
            io = nip1
            jo = 1
            je = nip1
          else if (iside .eq. 4) then
            i  = njp1
            io = nip1
            jo = 1
            je = nip1
          end if
          ioff = 1 + io*(i-1)

          do j = 1,je

            ij = ioff + jo*(j-1)
            dip(ij)  = swap_geo(pnti   )
            diu(ij)  = swap_geo(pnti+ 1)
            div(ij)  = swap_geo(pnti+ 2)
            dix(ij)  = swap_geo(pnti+ 3)
            rdip(ij) = swap_geo(pnti+ 4)
            rdiv(ij) = swap_geo(pnti+ 5)
            rdix(ij) = swap_geo(pnti+ 6)

            djp(ij)  = swap_geo(pnti+ 7)
            dju(ij)  = swap_geo(pnti+ 8)
            djv(ij)  = swap_geo(pnti+ 9)
            djx(ij)  = swap_geo(pnti+10)
            rdjp(ij) = swap_geo(pnti+11)
            rdju(ij) = swap_geo(pnti+12)
            rdjx(ij) = swap_geo(pnti+13)

            hijp(ij) = swap_geo(pnti+14)
            hjip(ij) = swap_geo(pnti+15)
            hkip(ij) = swap_geo(pnti+16)
            hkjp(ij) = swap_geo(pnti+17)

            hiju(ij) = swap_geo(pnti+18)
            hjiu(ij) = swap_geo(pnti+19)
            hkiu(ij) = swap_geo(pnti+20)
            hkju(ij) = swap_geo(pnti+21)

            hijv(ij) = swap_geo(pnti+22)
            hjiv(ij) = swap_geo(pnti+23)
            hkiv(ij) = swap_geo(pnti+24)
            hkjv(ij) = swap_geo(pnti+25)

            rp(ij) = swap_geo(pnti+26)
            ru(ij) = swap_geo(pnti+27)
            rv(ij) = swap_geo(pnti+28)

            didjp(ij)  = swap_geo(pnti+29)
            didju(ij)  = swap_geo(pnti+30)
            rdidju(ij) = swap_geo(pnti+31)
            djdiv(ij)  = swap_geo(pnti+32)
            rdjdiv(ij) = swap_geo(pnti+33)

            pnti = pnti + 34

          end do

        end if

      end if

      return
      end
