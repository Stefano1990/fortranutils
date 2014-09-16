
C=======================================================================
C CONTROL FOR SUCCESSIVE POINT OVERRELAXATION FOR A BLOCK
C=======================================================================
C To use this routine, set meth_smooth = 2
C Should then be used in place of sporlx on anything but
C ilev = nummgl (the coarsest grid).
C Necessary to check if ilev = nummgl in this routine also,
C sadly.

      subroutine sporlx2(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,niter,
     &  omega,soltol,work,res_norm,rhs_norm,err_norm,var_norm)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "les.inc"
#include "comm.inc"
#include "data_type.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"


C INPUT/OUTPUT PARAMETERS

      integer ilev
      integer niter
      real    f(*), b(*), residu(*)
      real    ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*)
      real    omega, soltol
      real    res_norm, rhs_norm, err_norm, var_norm
#ifdef DEBUG
      real    res_norm2(maxbl), rhs_norm2(maxbl)
      real    err_norm2(maxbl), var_norm2(maxbl)
#endif
      real    work(*)

C LOCAL

      integer iblk, iblk0,i,iline
      integer nip1, njp1, nkp1
      integer is, ie, js, je, ks, ke
      integer io, jo, ko
      integer p3
      integer isweep
      integer large_integer
      real    resnrm, errnrm, varnrm
      real    rnpts_all

      integer ierr
      real    red(2), red2(2)
      logical sporktest, swaptest


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C INITIAL ERROR

      resnrm = res_norm
      errnrm = err_norm
      varnrm = var_norm
      rnpts_all = 1./npts_all(ilev)

C SWEEP CONTROL

      isweep = 0
      do while (isweep.lt.niter .and. err_norm/var_norm.gt.soltol)

        isweep = isweep + 1

        do 10 iline = 1,7

          res_norm = 0.
          err_norm = 0.
          var_norm = 0.

C LOOP OVER BLOCKS
          iblk = 0
          do 15 iblk0 = 1,numbl_total
            if (bl2pe(iblk0) .eq. pe) then
              iblk = iblk + 1

              nip1 = mg_nijkp1(1,iblk,ilev)
              njp1 = mg_nijkp1(2,iblk,ilev)
              nkp1 = mg_nijkp1(3,iblk,ilev)

              io = 1
              jo = nip1
              ko = nip1*njp1

              is = nijkse(1,iblk,ilev)
              ie = nijkse(2,iblk,ilev)
              js = nijkse(3,iblk,ilev)
              je = nijkse(4,iblk,ilev)
              ks = nijkse(5,iblk,ilev)
              ke = nijkse(6,iblk,ilev)

              p3 = mg_off(iblk,ilev) + 1

              call spor(nip1,njp1,nkp1,io,jo,ko,is,ie,js,je,ks,ke,
     &             f(p3),b(p3),ap(p3),aw(p3),ae(p3),as(p3),an(p3),
     &             al(p3),ar(p3),omega,resnrm,errnrm,varnrm)

              res_norm = res_norm + resnrm
              err_norm = err_norm + errnrm
              var_norm = var_norm + varnrm

#ifdef DEBUG
              res_norm2(iblk) = resnrm
              err_norm2(iblk) = errnrm
              var_norm2(iblk) = varnrm
#endif

            endif
 15       continue

          call swap_sc(SV_P,ilev,f,work,1)

 10     end do


#ifdef DEBUG
        call usr_sum_double (err_norm2, red2(1), numbl_total,
     &  REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
        call usr_sum_double (var_norm2, red2(2), numbl_total,
     &  REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
        red(1) = err_norm
        red(2) = var_norm
        call MPI_ALLREDUCE(red, red2, 2, REAL_TYPE_MPI,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        err_norm = sqrt(red2(1)*rnpts_all)
        var_norm = sqrt(red2(2)*rnpts_all)
        if (var_norm .eq. 0.) var_norm = 1.

C END WHILE
      end do


#ifdef DEBUG
      call usr_sum_double (res_norm2, res_norm, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
      red(1) = res_norm
      call MPI_ALLREDUCE(red, res_norm, 1, REAL_TYPE_MPI,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      res_norm = sqrt(res_norm*rnpts_all)

      return
      end





