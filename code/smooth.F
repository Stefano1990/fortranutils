C=======================================================================
C SMOOTHING CONTROL
C=======================================================================

      subroutine smooth(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,nsmooth,
     &  omega,soltol,work,offwork,resnrm,rhsnrm,errnrm,varnrm)

      implicit none

#include "dim.inc"
#include "comm.inc"

#include "cb.control.inc"
#include "cb.comm.inc"

      integer ilev
      integer nsmooth
      real    ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*)
      real    f(*), b(*), residu(*)
      real    omega, soltol
      real    resnrm, rhsnrm, errnrm, varnrm
      integer offwork(*)
      real    work(*)

C LOCAL

      integer ismooth

      real    rdum(3)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      if (ilev.eq.1 .or. ilev.eq.nummgl) then
        call norms(ilev,1,ap,aw,ae,as,an,al,ar,f,b,residu,
     &       resnrm,rhsnrm,errnrm,varnrm)
      else
        errnrm = 1.
        varnrm = 1.
      end if

      ismooth = 0

C DEBUG INFORMATION OUTPUT

c     call global_sum(ilev,1,b,rdum(1),rdum(2),rdum(3))
c     write(*,*) '************ right hand sum: ',rdum(1)

c     if (master) then
c       write(*,'(''*** Info: Begin in smooth: ilev='',i2,
c    &            '' ||res||='',e10.4,)') ilev,resnrm
c     end if

 10   if (errnrm/varnrm.gt.soltol .and. ismooth.lt.nsmooth) then

        ismooth = ismooth + 1

C ALTERNATING LINE

C SUCCESSIVE POINT OVERRELAXATION

        if (meth_smooth .eq. 1) then

          call sporlx(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,1,omega,
     &         soltol,work,offwork,resnrm,rhsnrm,errnrm,varnrm)

        else if (meth_smooth .eq. 2) then

c         call rbsporlx(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,1,omega,
c    &         soltol,work,offwork,resnrm,rhsnrm,errnrm,varnrm)

          call sporlx2(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,1,omega,
     &         soltol,work,resnrm,rhsnrm,errnrm,varnrm)

        else

          write(*,'(''*** Error: Smoothing method No. '',i4,
     &              '' not implemented yet'')') meth_smooth
          stop

        end if

 15     continue

C CALCULATE THE NORMS AND RESIDUAL AFTER SMOOTHING

        if (ismooth .eq. nsmooth) then
          call norms(ilev,1,ap,aw,ae,as,an,al,ar,f,b,residu,
     &         resnrm,rhsnrm,errnrm,varnrm)

C DEBUG INFORMATION OUTPUT
c         if (master) then
c           write(*,'(''*** Info: End   in smooth: ilev='',i2,
c    &                '' ||res||='',e10.4,'' nsmooth='',i4)')
c    &                ilev,resnrm,ismooth
c         end if

        end if

        goto 10

      end if

      return
      end
