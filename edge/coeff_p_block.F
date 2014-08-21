C=======================================================================
C COEFFICIENTS FOR PRESSURE EQUATION
C=======================================================================

      subroutine coeff_p_block(nip1,njp1,nkp1,jo,ko,
     &                   diu, djv, rdiv, rdju, didjp, rp,
     &                   dkp, dkw,
     &                   ap,aw,ae,as,an,al,ar,bc)

      implicit none

#include "bc_type.inc"

C ARGUMENTS

      integer      nip1, njp1, nkp1
      integer      jo, ko

      real         diu(nip1,njp1)
      real         djv(nip1,njp1)
      real         rdiv(nip1,njp1)
      real         rdju(nip1,njp1)
      real         didjp(nip1,njp1)
      real         rp(nip1,njp1)

      real         dkp(nkp1)
      real         dkw(nkp1)

      real         ap(nip1*njp1*nkp1)
      real         aw(nip1*njp1*nkp1)
      real         ae(nip1*njp1*nkp1)
      real         as(nip1*njp1*nkp1)
      real         an(nip1*njp1*nkp1)
      real         al(nip1*njp1*nkp1)
      real         ar(nip1*njp1*nkp1)

      integer      bc(6)

C LOCAL

      integer      i, j, k
      integer      ni, nj, nk
      integer      ioff, joff, koff, ijoff
      integer      ijk
      integer      ijkw, ijks, ijkl
      real         cef, cnf, crf


C COMMON

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------

      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C-----------------------------------------------------------------------
C DIAGNOSTIC HEADING
C-----------------------------------------------------------------------
 
C-----------------------------------------------------------------------
C SET COEFFICIENTS
C-----------------------------------------------------------------------

C INITIALIZE THE STARTING SURFACES

C WEST

      ioff = 1
      do k = 2,nk
        koff = ko*(k-1) + ioff
        do j = 2,nj
          ijk = jo*(j-1) + koff
          ae(ijk) = 0.
        end do
      end do

C SOUTH

      joff = 0
      do k = 2,nk
        koff = ko*(k-1) + joff
        do i = 2,ni
          ijk = i + koff
          an(ijk) = 0.
        end do
      end do

C LEFT

      koff = 0
      do j = 2,nj
        joff = jo*(j-1) + koff
        do i = 2,ni
          ijk = i + joff
          ar(ijk) = 0.
        end do
      end do

C LOOP THROUGH IJ SURFACE

      do j=2,nj
        do i=2,ni

          cef = rdju(i+1,j) / diu(i+1,j)
          cnf = rdiv(i,j+1) / djv(i,j+1)
          crf = didjp(i,j)  / rp(i,j)

C LOOP THROUGH K   ?????? SHOULD HAVE A STRIDE OF 1! ???????

          ijoff = i + jo*(j-1) - ko

          do k=2,nk

            ijk  = ijoff + ko*k
            ijkw = ijk - 1
            ijks = ijk - jo
            ijkl = ijk - ko

            aw(ijk) = ae(ijkw)
            ae(ijk) = cef*dkp(k)

            as(ijk) = an(ijks)
            an(ijk) = cnf*dkp(k)

            al(ijk) = ar(ijkl)
            ar(ijk) = crf/dkw(k+1)

            ap(ijk) = aw(ijk) + ae(ijk) + as(ijk) + an(ijk) +
     &                al(ijk) + ar(ijk)

          end do

        end do
      end do

C SET UP THE STARTING/ENDING SURFACES. THE COEFFICIENT IS SET TO ZERO
C IF THE CORRESPONDING B.C. IS A WALL CONDITION (THIS IMPLIES NEUMANN
C CONDITION FOR THE PRESSURE)

C WEST

      i    = 2
      ioff = 2
      if (bc(1).eq.BCT_PERIODIC .or. bc(1).eq.BCT_INTERFACE) then
        do k = 2,nk
          koff = ko*(k-1) + ioff
          do j = 2,nj
            ijk = jo*(j-1) + koff
            aw(ijk) = rdju(i,j)/diu(i,j)*dkp(k)
            ap(ijk) = ap(ijk) + aw(ijk)
          end do
        end do
      end if

C EAST

      ioff = ni
      if (bc(2).ne.BCT_PERIODIC .and. bc(2).ne.BCT_INTERFACE) then
        do k = 2,nk
          koff = ko*(k-1) + ioff
          do j = 2,nj
            ijk = jo*(j-1) + koff
            ap(ijk) = ap(ijk) - ae(ijk)
            ae(ijk) = 0.
          end do
        end do
      end if

C SOUTH

      j    = 2
      joff = jo
      if (bc(3).eq.BCT_PERIODIC .or. bc(3).eq.BCT_INTERFACE) then
        do k = 2,nk
          koff = ko*(k-1) + joff
          do i = 2,ni
            ijk = i + koff
            as(ijk) = rdiv(i,j)/djv(i,j)*dkp(k)
            ap(ijk) = ap(ijk) + as(ijk)
          end do
        end do
      end if

C NORTH

      joff = jo*(nj-1)
      if (bc(4).ne.BCT_PERIODIC .and. bc(4).ne.BCT_INTERFACE) then
        do k = 2,nk
          koff = ko*(k-1) + joff
          do i = 2,ni
            ijk = i + koff
            ap(ijk) = ap(ijk) - an(ijk)
            an(ijk) = 0.
          end do
        end do
      end if

C LEFT

      k    = 2
      koff = ko
      if (bc(5).eq.BCT_PERIODIC .or. bc(5).eq.BCT_INTERFACE) then
        do j = 2,nj
          joff = jo*(j-1) + koff
          do i = 2,ni
            ijk = i + joff
            al(ijk) = didjp(i,j)/(rp(i,j)*dkw(k))
            ap(ijk) = ap(ijk) + al(ijk)
          end do
        end do
      end if

C RIGHT

      koff = ko*(nk-1)
      if (bc(6).ne.BCT_PERIODIC .and. bc(6).ne.BCT_INTERFACE) then
        do j = 2,nj
          joff = jo*(j-1) + koff
          do i = 2,ni
            ijk = i + joff
            ap(ijk) = ap(ijk) - ar(ijk)
            ar(ijk) = 0.
          end do
        end do
      end if

      return
      end
      
