C=======================================================================
C CONTROL FOR SUCCESSIVE POINT OVERRELAXATION FOR A BLOCK
C=======================================================================

      subroutine sporlx(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,niter,
     &  omega,soltol,work,offwork,res_norm,rhs_norm,err_norm,var_norm)

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
      integer offwork(*)
#ifdef DEBUG
      real    res_norm2(maxbl), rhs_norm2(maxbl)
      real    err_norm2(maxbl), var_norm2(maxbl)
#endif
      real    work(*)

C LOCAL

      integer iblk
      integer d1, nip1, njp1, nkp1
      integer is, ie, js, je, ks, ke
      integer io, jo, ko
      integer p3
      integer isweep, iline, jumble, jumst, d(3), jumfor
      integer large_integer
      real    resnrm, errnrm, varnrm
      real    rnpts_all

      integer ierr
      real    red(2),red2(2)
      integer status(MPI_STATUS_SIZE), i



      data    large_integer/99999/, jumble/0/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C INITIAL ERROR

      resnrm = res_norm
      errnrm = err_norm
      varnrm = var_norm
      rnpts_all = 1./npts_all(ilev)


C SET SWEEP DIRECTION (THE DIRECTION OF A SWEEP IS THE DIRECTION OF
C THE COORDINATE IN WHICH THE POINT SWEEPING MOVES)

      if(mod(jumble,16).ge.8) then
        jumfor = 1
      else
        jumfor = 2
      end if

      jumst = mod(jumble,24)/8 + 1

      d(1) = jumst
      d(2) = cyclic(mod(jumfor+1,2)+1,jumst)
      d(3) = cyclic(mod(jumfor  ,2)+1,jumst)

C SWEEP CONTROL

      isweep = 0
      do while (isweep.lt.niter .and. err_norm/var_norm.gt.soltol)

        isweep = isweep + 1

C LOOP OVER DIRECTIONS

        do 10 iline = 1,3

          d1 = d(iline)
c         npts_all = 0.
          res_norm = 0.
          err_norm = 0.
          var_norm = 0.

C LOOP OVER BLOCKS

          do 15 iblk = 1,numbl

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

            if (d1 .eq. 1) then
              call spor(nip1,njp1,nkp1,io,jo,ko,is,ie,js,je,ks,ke,
     &             f(p3),b(p3),ap(p3),aw(p3),ae(p3),as(p3),an(p3),
     &             al(p3),ar(p3),omega,resnrm,errnrm,varnrm)
            else if (d1 .eq. 2) then
              call spor(nip1,njp1,nkp1,jo,io,ko,js,je,is,ie,ks,ke,
     &             f(p3),b(p3),ap(p3),as(p3),an(p3),aw(p3),ae(p3),
     &             al(p3),ar(p3),omega,resnrm,errnrm,varnrm)
            else
              call spor(nip1,njp1,nkp1,ko,io,jo,ks,ke,is,ie,js,je,
     &             f(p3),b(p3),ap(p3),al(p3),ar(p3),aw(p3),ae(p3),
     &             as(p3),an(p3),omega,resnrm,errnrm,varnrm)
            end if

            res_norm = res_norm + resnrm
            err_norm = err_norm + errnrm
            var_norm = var_norm + varnrm

#ifdef DEBUG
            res_norm2(iblk) = resnrm
            err_norm2(iblk) = errnrm
            var_norm2(iblk) = varnrm
#endif

 15       continue

C APPLICATION OF PERIODIC B.C. AND EXCHANGE INFORMATION AMONG INTERFACE
C POINTS
          call swap_sc(SV_P,ilev,f,work,1)

C!!!! debug

c         write(*,'(''iter, norms: '',i4,2x,4(1x,e15.8))')
c    &    isweep,res_norm,rhs_norm,err_norm,var_norm

 10     continue

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
c       if(master)write(6,*)'Error ',err_norm

C END WHILE

      end do

#ifdef DEBUG
      call usr_sum_double (res_norm2, res_norm, numbl_total,
     &REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
      red = res_norm
      call MPI_ALLREDUCE(red, res_norm, 1, REAL_TYPE_MPI,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      res_norm = sqrt(res_norm*rnpts_all)

      jumble = jumble + 13

c     SET jumble BACK TO A SMALL VALUE WHEN IT'S TOO LARGE

      if (jumble .gt. large_integer) then
        jumble = mod(jumble,large_integer)
      end if

      return
      end

C=======================================================================
C SUCCESSIVE POINT OVERRELAXATION FOR A BLOCK
C=======================================================================

C 

      subroutine spor(nip1,njp1,nkp1,io,jo,ko,is,ie,js,je,ks,ke,
     &  f,src,ap,aw,ae,as,an,al,ar,omega,resnrm,errnrm,varnrm)

C INPUT/OUTPUT 

      implicit none

      integer is, ie, js, je, ks, ke
      integer nip1, njp1, nkp1
      integer io, jo, ko, koff, jkoff
      real    ap(nip1*njp1*nkp1), aw(nip1*njp1*nkp1),
     &        ae(nip1*njp1*nkp1), as(nip1*njp1*nkp1),
     &        an(nip1*njp1*nkp1), al(nip1*njp1*nkp1),
     &        ar(nip1*njp1*nkp1), f(nip1*njp1*nkp1),
     &        src(nip1*njp1*nkp1)
      real    omega, resnrm, errnrm, varnrm

C LOCAL

      integer i, j, k, ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr
      real    dx, f_old, res

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      resnrm = 0.
      errnrm = 0.
      varnrm = 0.

C LOOP OVER GRID POINTS TO CARRY OUT GAUSS-SEIDEL SUCCESSIVE POINT
C OVERRELAXATION (SPOR)

      do 10 k = ks,ke
        koff = 1 + ko*(k-1)
        do 15 j = js,je
          jkoff = jo*(j-1) + koff
          do 20 i = is,ie

            ijk  = io*(i-1) + jkoff
            ijkw = ijk - io
            ijke = ijk + io
            ijks = ijk      - jo
            ijkn = ijk      + jo
            ijkl = ijk           - ko
            ijkr = ijk           + ko

            f_old  = f(ijk)
            dx = (src(ijk) +
     &           aw(ijk)*f(ijkw) + ae(ijk)*f(ijke) +
     &           as(ijk)*f(ijks) + an(ijk)*f(ijkn) +
     &           al(ijk)*f(ijkl) + ar(ijk)*f(ijkr))/
     &           ap(ijk) - f_old
            f(ijk) = f_old + omega*dx

            dx     = f(ijk) - f_old
            res    = dx*ap(ijk)
            resnrm = resnrm + res*res
            errnrm = errnrm + dx*dx
            varnrm = varnrm + f(ijk)*f(ijk)

 20       continue
 15     continue
 10   continue

      return
      end
