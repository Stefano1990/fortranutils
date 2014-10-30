C=======================================================================
C CALCULATE NORMS FOR A BLOCK
C=======================================================================

      subroutine norms_block(ap,aw,ae,as,an,al,ar,x,b,residu,nip1,njp1,
     &  nkp1,nis,nie,njs,nje,nks,nke,res_norm,rhs_norm,err_norm,
     &  var_norm)

      implicit none

C INPUT/OUTPUT PARAMETERS

      integer nip1, njp1, nkp1
      integer nis, nie, njs, nje, nks, nke
      real    ap(nip1*njp1*nkp1), aw(nip1*njp1*nkp1),
     &        ae(nip1*njp1*nkp1), as(nip1*njp1*nkp1),
     &        an(nip1*njp1*nkp1), al(nip1*njp1*nkp1),
     &        ar(nip1*njp1*nkp1), x(nip1*njp1*nkp1),
     &        b(nip1*njp1*nkp1), residu(nip1*njp1*nkp1)
      real    res_norm, rhs_norm, err_norm, var_norm

C LOCAL QUANTITIES

      integer i, j, k, ijk
      integer io, jo, ko, koff, jkoff
      real    res, rhs, err, var

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      res_norm = 0.
      rhs_norm = 0.
      err_norm = 0.
      var_norm = 0.

      io   = 1
      jo   = nip1
      ko   = nip1*njp1

      do k = nks,nke
        koff = ko*(k-1)
        do j = njs,nje
          jkoff = jo*(j-1) + koff
          do i = nis,nie
            ijk = i + jkoff

            rhs = b(ijk)
            var = x(ijk)
            res = b(ijk) - ap(ijk)*var +
     &            aw(ijk)*x(ijk-io) + ae(ijk)*x(ijk+io) +
     &            as(ijk)*x(ijk-jo) + an(ijk)*x(ijk+jo) +
     &            al(ijk)*x(ijk-ko) + ar(ijk)*x(ijk+ko)
            err = res/ap(ijk)

            res_norm = res_norm + res*res
            rhs_norm = rhs_norm + rhs*rhs
            err_norm = err_norm + err*err
            var_norm = var_norm + var*var

            residu(ijk) = res

          end do
        end do
      end do

      return
      end
