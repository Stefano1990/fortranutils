C=======================================================================
C RESTRICTION FOR A BLOCK
C 1. - imult, jmult, kmult ARE ASSUMED TO BE 2 AT THIS MOMENT
C=======================================================================

      subroutine rest_block(ff,fc,
     &  nip1c,njp1c,nkp1c,joc,koc,
     &  nip1f,njp1f,nkp1f,jof,kof,
     &  imult,jmult,kmult)

      implicit none

      integer nip1c, njp1c, nkp1c, joc, koc
      integer nip1f, njp1f, nkp1f, jof, kof
      integer imult, jmult, kmult
      real    fc(nip1c*njp1c*nkp1c), ff(nip1f*njp1f*nkp1f)

C LOCAL

      integer i, j, k
      integer ijkc, ijkf
      integer koffc, jkoffc, kofff, jkofff
      integer offm11, off1m1, offmm1, off11m, offm1m, off1mm, offmmm
      integer iinc, jinc, kinc
      real    fac

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      fac = 1.0

      if (imult.eq.1) then
        iinc  = 0
        fac   = 0.5*fac
      else
        iinc  = -1
      end if

      if (jmult.eq.1) then
        jinc  = 0
        fac   = 0.5*fac
      else
        jinc  = -jof
      end if

      if (kmult.eq.1) then
        kinc  = 0
        fac   = 0.5*fac
      else
        kinc  = -kof
      end if

      offm11 = iinc
      off1m1 =        jinc
      offmm1 = iinc + jinc
      off11m =               kinc
      offm1m = iinc        + kinc
      off1mm =        jinc + kinc
      offmmm = iinc + jinc + kinc

      do 10 k = 2,nkp1c-1
        koffc = koc*(k-1)
        kofff = kof*(kmult*k-2)
        do 15 j = 2,njp1c-1
          jkoffc = joc*(j-1)       + koffc
          jkofff = jof*(jmult*j-2) + kofff
          do 20 i = 2,nip1c-1
            ijkc = i         + jkoffc
            ijkf = imult*i-1 + jkofff
            fc(ijkc) = fac*(ff(ijkf       ) + ff(ijkf+offm11) +
     &                      ff(ijkf+off1m1) + ff(ijkf+offmm1) +
     &                      ff(ijkf+off11m) + ff(ijkf+offm1m) +
     &                      ff(ijkf+off1mm) + ff(ijkf+offmmm))
 20       continue
 15     continue
 10   continue
      
      return
      end
