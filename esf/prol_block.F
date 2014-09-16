C=======================================================================
C PROLONGATION
C - 1. TRILINEAR INTERPOLATION
C - 2. MODIFICATION TO PROLONGATION WITH NEUMANN BOUNDARY CONDITION
C - 3. IT IS ASSUMED THAT THE COARSE GRID COVERS THE FINE GRID
C - 4. THE CORRECTION IS CALCULATED FOR THE VIRTUAL FINE GRID
C - 5. THIS CODE WORKS FOR CASES WITHOUT SPACE DIMENSION DEGENERATION
C=======================================================================

      subroutine prol_block(ff,fc,cr,
     &  nip1c,njp1c,nkp1c,nip1f,njp1f,nkp1f,
     &  imult,jmult,kmult,bc)

      implicit none

#include "dim.inc"
#include "bc_type.inc"

#include "cb.lgrid.inc"
#include "cb.lcontrol.inc"

      integer nip1c , njp1c , nkp1c , nip1f, njp1f, nkp1f
      real    fc(nip1c*njp1c*nkp1c), cr(nip1f*njp1f*nkp1f),
     &        ff(nip1f*njp1f*nkp1f)
      integer imult, jmult, kmult
      integer bc(6)

C-----------------------------------------------------------------------

      integer i_c, j_c, k_c, i_f, j_f, k_f
      integer ijkc, ijkf
      integer joc, koc, jkoffc, koffc
      integer jof, kof, jkofff, kofff
      integer ioc_tmp, joc_tmp, koc_tmp
      integer eoffc, noffc, roffc, neoffc, nroffc, eroffc, enroffc,
     &        eofff, nofff, rofff, neofff, nrofff, erofff, enrofff
      integer inext, ioffc, iside
      integer ijkoc(3) , nip1c_tmp, njp1c_tmp, nkp1c_tmp, nijkp1c(3)
      integer ijkof(3)
      integer d1, d2, d3
      integer is, ie, js, je, ks, ke


      real    c0
      real    c1i, c2i, c3i
      real    c1s, c2s, c3s
      real    c1e, c2e

#ifdef T3E
      data    c0/0.015625/
      data    c1i, c2i, c3i/0.421875, 0.140625 , 0.046875/
      data    c1s, c2s, c3s/0.5625, 0.1875, 0.0625/
      data    c1e, c2e/0.75, 0.25/
#else
      data    c0/0.015625d0/
      data    c1i, c2i, c3i/0.421875d0, 0.140625d0 , 0.046875d0/
      data    c1s, c2s, c3s/0.5625d0, 0.1875d0, 0.0625d0/
      data    c1e, c2e/0.75d0, 0.25d0/
#endif

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      joc         = nip1c
      koc         = nip1c*njp1c
      ijkoc(1)    = 1
      ijkoc(2)    = joc
      ijkoc(3)    = koc
      nijkp1c(1)  = nip1c
      nijkp1c(2)  = njp1c
      nijkp1c(3)  = nkp1c

      jof         = nip1f
      kof         = nip1f*njp1f
      ijkof(1)    = 1
      ijkof(2)    = jof
      ijkof(3)    = kof

      eoffc   = 1
      noffc   =     joc
      roffc   =           koc 
      neoffc  = 1 + joc
      nroffc  =     joc + koc
      eroffc  = 1       + koc
      enroffc = 1 + joc + koc

      eofff   = 1
      nofff   =     jof
      rofff   =           kof 
      neofff  = 1 + jof
      nrofff  =     jof + kof
      erofff  = 1       + kof
      enrofff = 1 + jof + kof

C-----------------------------------------------------------------------
C COPY THE CORRECTION VALUES ON THE POINTS NEXT TO THE BOUNDARY TO
C THE POINTS ON THE BOUNDARY, IF THE BOUNDARY IS NOT AN INTERFACE.
C THIS IS TO SIMPLIFY THE CODING FOR THE PROLONGATION.
C-----------------------------------------------------------------------

      do 10 iside = 1,6

        d1 = (iside+1)/2
        d2 = cyclic1(1,d1)
        d3 = cyclic1(2,d1)

        nip1c_tmp = nijkp1c(d1)
        njp1c_tmp = nijkp1c(d2)
        nkp1c_tmp = nijkp1c(d3)
        ioc_tmp = ijkoc(d1)
        joc_tmp = ijkoc(d2)
        koc_tmp = ijkoc(d3)

        if (bc(iside) .ne. BCT_INTERFACE) THEN

          if (mod(iside,2).eq.1) then
            ioffc = 1
            inext =  ioc_tmp
          else
            ioffc = 1 + ioc_tmp*(nip1c_tmp-1)
            inext = -ioc_tmp
          end if
          do 15 k_c = 1,nkp1c_tmp
            koffc = ioffc + koc_tmp*(k_c-1)
            do 20 j_c = 1,njp1c_tmp
              ijkc = joc_tmp*(j_c-1) + koffc
              fc(ijkc) = fc(ijkc+inext)
 20         continue
 15       continue

        end if

 10   continue

C-----------------------------------------------------------------------
C PROLONGATE COARSE GRID CORRECTION
C-----------------------------------------------------------------------

      do i_c = 1,nip1f*njp1f*nkp1f
        cr(i_c) = 0.
      end do

      do 25 k_c = 1,nkp1c-1
        k_f   = kmult*k_c - 1
        koffc = koc*(k_c-1)
        kofff = kof*(k_f-1)
        do 30 j_c = 1,njp1c-1
          j_f    = jmult*j_c - 1
          jkoffc = joc*(j_c-1) + koffc
          jkofff = jof*(j_f-1) + kofff
          do 35 i_c = 1,nip1c-1
            i_f  = imult*i_c - 1
            ijkc = i_c + jkoffc
            ijkf = i_f + jkofff

            cr(ijkf        ) = cr(ijkf        ) +
     &        c1i*fc(ijkc       ) + c2i*fc(ijkc+eoffc  ) +
     &        c2i*fc(ijkc+noffc ) + c2i*fc(ijkc+roffc  ) +
     &        c3i*fc(ijkc+neoffc) + c3i*fc(ijkc+nroffc ) +
     &        c3i*fc(ijkc+eroffc) + c0 *fc(ijkc+enroffc)
            cr(ijkf+eofff  ) = cr(ijkf+eofff  ) +
     &        c2i*fc(ijkc       ) + c1i*fc(ijkc+eoffc  ) +
     &        c3i*fc(ijkc+noffc ) + c3i*fc(ijkc+roffc  ) +
     &        c2i*fc(ijkc+neoffc) + c0 *fc(ijkc+nroffc ) +
     &        c2i*fc(ijkc+eroffc) + c3i*fc(ijkc+enroffc)
            cr(ijkf+nofff  ) = cr(ijkf+nofff  ) +
     &        c2i*fc(ijkc       ) + c3i*fc(ijkc+eoffc  ) +
     &        c1i*fc(ijkc+noffc ) + c3i*fc(ijkc+roffc  ) +
     &        c2i*fc(ijkc+neoffc) + c2i*fc(ijkc+nroffc ) +
     &        c0 *fc(ijkc+eroffc) + c3i*fc(ijkc+enroffc)
            cr(ijkf+rofff  ) = cr(ijkf+rofff  ) +
     &        c2i*fc(ijkc       ) + c3i*fc(ijkc+eoffc  ) +
     &        c3i*fc(ijkc+noffc ) + c1i*fc(ijkc+roffc  ) +
     &        c0 *fc(ijkc+neoffc) + c2i*fc(ijkc+nroffc ) +
     &        c2i*fc(ijkc+eroffc) + c3i*fc(ijkc+enroffc)
            cr(ijkf+neofff ) = cr(ijkf+neofff ) +
     &        c3i*fc(ijkc       ) + c2i*fc(ijkc+eoffc  ) +
     &        c2i*fc(ijkc+noffc ) + c0 *fc(ijkc+roffc  ) +
     &        c1i*fc(ijkc+neoffc) + c3i*fc(ijkc+nroffc ) +
     &        c3i*fc(ijkc+eroffc) + c2i*fc(ijkc+enroffc)
            cr(ijkf+nrofff ) = cr(ijkf+nrofff ) +
     &        c3i*fc(ijkc       ) + c0 *fc(ijkc+eoffc  ) +
     &        c2i*fc(ijkc+noffc ) + c2i*fc(ijkc+roffc  ) +
     &        c3i*fc(ijkc+neoffc) + c1i*fc(ijkc+nroffc ) +
     &        c3i*fc(ijkc+eroffc) + c2i*fc(ijkc+enroffc)
            cr(ijkf+erofff ) = cr(ijkf+erofff ) +
     &        c3i*fc(ijkc       ) + c2i*fc(ijkc+eoffc  ) +
     &        c0 *fc(ijkc+noffc ) + c2i*fc(ijkc+roffc  ) +
     &        c3i*fc(ijkc+neoffc) + c3i*fc(ijkc+nroffc ) +
     &        c1i*fc(ijkc+eroffc) + c2i*fc(ijkc+enroffc)
            cr(ijkf+enrofff) = cr(ijkf+enrofff) +
     &        c0 *fc(ijkc       ) + c3i*fc(ijkc+eoffc  ) +
     &        c3i*fc(ijkc+noffc ) + c3i*fc(ijkc+roffc  ) +
     &        c2i*fc(ijkc+neoffc) + c2i*fc(ijkc+nroffc ) +
     &        c2i*fc(ijkc+eroffc) + c1i*fc(ijkc+enroffc)

 35       continue
 30     continue
 25   continue

C-----------------------------------------------------------------------
C CARRY OUT COARSE GRID CORRECTION
C-----------------------------------------------------------------------

      if (bc(1) .eq. BCT_INTERFACE) then
        is = 1
      else
        is = 2
      end if
      if (bc(2) .eq. BCT_INTERFACE) then
        ie = nip1f
      else
        ie = nip1f - 1
      end if
      if (bc(3) .eq. BCT_INTERFACE) then
        js = 1
      else
        js = 2
      end if
      if (bc(4) .eq. BCT_INTERFACE) then
        je = njp1f
      else
        je = njp1f - 1
      end if
      if (bc(5) .eq. BCT_INTERFACE) then
        ks = 1
      else
        ks = 2
      end if
      if (bc(6) .eq. BCT_INTERFACE) then
        ke = nkp1f
      else
        ke = nkp1f - 1
      end if

      do 40 k_f = ks,ke
        kofff = kof*(k_f-1)
        do 45 j_f = js,je
          jkofff = jof*(j_f-1) + kofff
          do 50 i_f = is,ie
            ijkf = i_f + jkofff
            ff(ijkf) = ff(ijkf) + cr(ijkf)
 50       continue
 45     continue
 40   continue

      return
      end
