C=======================================================================
C CALCULATE NORMS FOR ALL BLOCKS
C=======================================================================

      subroutine norms(ilev,mg_or_sg,ap,aw,ae,as,an,al,ar,x,b,residu,
     &  res_norm,rhs_norm,err_norm,var_norm)

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"

C COMMON BLOCKS
#include "cb.lgrid.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"


C INPUT/OUTPUT PARAMETERS

      integer ilev, mg_or_sg
      real    ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*),
     &        x(*), b(*), residu(*)
      real    res_norm, rhs_norm, err_norm, var_norm

C LOCAL QUANTITIES

      integer iblk
      integer nip1, njp1, nkp1, is, ie, js, je, ks, ke
      integer p3, p3x
      real    resnrm, rhsnrm, errnrm, varnrm
      real    red1(4), red2(4), rnpts_all
      integer ierr

      integer io,jo,ko,ijk,i,j,k,koff,jkoff

      integer status(MPI_STATUS_SIZE)
#ifdef DEBUG
      real    res_norm2(maxbl), rhs_norm2(maxbl)
      real    err_norm2(maxbl), var_norm2(maxbl)
#endif

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      res_norm = 0.
      rhs_norm = 0.
      err_norm = 0.
      var_norm = 0.


      do 10 iblk = 1,numbl

        p3 = mg_off(iblk,ilev) + 1
        if (mg_or_sg .eq. 1) then
          p3x = p3
        else
          p3x = sg_off(iblk) + 1
        end if

        nip1 = mg_nijkp1(1,iblk,ilev)
        njp1 = mg_nijkp1(2,iblk,ilev)
        nkp1 = mg_nijkp1(3,iblk,ilev)

        is = nijkse(1,iblk,ilev)
        ie = nijkse(2,iblk,ilev)
        js = nijkse(3,iblk,ilev)
        je = nijkse(4,iblk,ilev)
        ks = nijkse(5,iblk,ilev)
        ke = nijkse(6,iblk,ilev)

        call norms_block(ap(p3),aw(p3),ae(p3),as(p3),an(p3),al(p3),
     &       ar(p3),x(p3x),b(p3),residu(p3),nip1,njp1,nkp1,
     &       is,ie,js,je,ks,ke,resnrm,rhsnrm,errnrm,varnrm)

        res_norm = res_norm + resnrm
        rhs_norm = rhs_norm + rhsnrm
        err_norm = err_norm + errnrm
        var_norm = var_norm + varnrm

#ifdef DEBUG
        res_norm2(iblk) = resnrm
        rhs_norm2(iblk) = rhsnrm
        err_norm2(iblk) = errnrm
        var_norm2(iblk) = varnrm
#endif

 10   continue

#ifdef DEBUG
      call usr_sum_double (res_norm2, res_norm, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
      call usr_sum_double (rhs_norm2, rhs_norm, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
      call usr_sum_double (err_norm2, err_norm, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
      call usr_sum_double (var_norm2, var_norm, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
      red1(1) = res_norm
      red1(2) = rhs_norm
      red1(3) = err_norm
      red1(4) = var_norm
      call MPI_ALLREDUCE(red1, red2, 4, REAL_TYPE_MPI,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
      res_norm = red2(1)
      rhs_norm = red2(2)
      err_norm = red2(3)
      var_norm = red2(4)
#endif

      rnpts_all = 1./npts_all(ilev)
      res_norm = sqrt(res_norm*rnpts_all)
      rhs_norm = sqrt(rhs_norm*rnpts_all)
      err_norm = sqrt(err_norm*rnpts_all)
      var_norm = sqrt(var_norm*rnpts_all)
      if (var_norm .eq. 0.) var_norm = 1.


      return
      end

