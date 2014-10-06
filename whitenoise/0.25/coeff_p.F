C=======================================================================
C LOOP THROUGH BLOCKS TO CALCULATE COEFFICIENTS FOR PRESSURE EQUATION
C=======================================================================

      subroutine coeff_p(diu, djv, rdiv, rdju, didjp, rp,
     &                   dkp, dkw,
     &                   ap,aw,ae,as,an,al,ar)

      implicit none

C PARAMETERS & INCLUDE FILES

#include "dim.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"

C ARGUMENTS

      real         diu(*), djv(*), rdiv(*), rdju(*), didjp(*)
      real         rp(*)

      real         dkp(*), dkw(*)

      real         ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*)

C LOCAL

      integer      ilev, iblk
      integer      nip1, njp1, nkp1
      integer      jo, ko
      integer      p1, p2, p3
      real         rdum(2)


C COMMON

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do 10 ilev = 1,nummgl
        do 10 iblk = 1,numbl

          nip1 = mg_nijkp1(1,iblk,ilev)
          njp1 = mg_nijkp1(2,iblk,ilev)
          nkp1 = mg_nijkp1(3,iblk,ilev)

          jo = nip1
          ko = nip1*njp1

C PINTERS

          p1 = mg_off_1d_k(iblk,ilev)  + 1
          p2 = mg_off_2d_ij(iblk,ilev) + 1
          p3 = mg_off(iblk,ilev)       + 1

          call coeff_p_block(nip1,njp1,nkp1,jo,ko,
     &         diu(p2),djv(p2),rdiv(p2),rdju(p2),didjp(p2),rp(p2),
     &         dkp(p1),dkw(p1),
     &         ap(p3),aw(p3),ae(p3),as(p3),an(p3),al(p3),ar(p3),
     &         bl_bc(1,iblk))

 15     continue

C!!!! debug

c       call file_3(20,'ap',1,ap,ilev,1,rdum(1),rdum(2))
c       call file_3(20,'aw',1,aw,ilev,1,rdum(1),rdum(2))
c       call file_3(20,'ae',1,ae,ilev,1,rdum(1),rdum(2))
c       call file_3(20,'as',1,as,ilev,1,rdum(1),rdum(2))
c       call file_3(20,'an',1,an,ilev,1,rdum(1),rdum(2))
c       call file_3(20,'al',1,al,ilev,1,rdum(1),rdum(2))
c       call file_3(20,'ar',1,ar,ilev,1,rdum(1),rdum(2))

10    continue

      return
      end
