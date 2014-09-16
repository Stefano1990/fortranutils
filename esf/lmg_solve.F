C=======================================================================
C  LINEAR MULTIGRID CONTROL
C=======================================================================

      subroutine lmg_solve(ap,aw,ae,as,an,al,ar,f,b,residu,ia,ra)

C-----------------------------------------------------------------------

      implicit none

#include "mpif.h"

#include "dim.inc"
#include "comm.inc"
#include "const.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.address.inc"
#include "cb.comm.inc"

c   Input/Output Parameters
 
      integer   ia(*)
      real      ra(*)
      real      ap(*), aw(*), ae(*), as(*), an(*), al(*), ar(*)
      real      f(*), b(*), residu(*)

c   Local Quantities
 
      integer   ilev, ierr
      integer   mg_gamma, maxgam
      parameter (maxgam=10)
      integer   nummg_l(maxmgl)
      real      rho_last, norm_last
      logical   fcycle, vcycle, goto_fine_grid
 
C-----------------------------------------------------------------------
 
C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------

      vcycle = meth_cycle.eq.1
      fcycle = meth_cycle.eq.2
 
      if(vcycle) then
        mg_gamma = 1
      else
        mg_gamma = 2
      end if
      tolc  = 0.
      if(nummgl.eq.1) tolc = tolf

C-----------------------------------------------------------------------
C COMPUTE THE ACCURACY TOLERANCE FOR THE SOLUTION ON THE FINEST GRID
C AND SET INITIAL CONTROL QUANTITIES
C-----------------------------------------------------------------------
 
      call norms(1,1,ap,aw,ae,as,an,al,ar,f,b,residu,
     &     rnorm(1),fnorm(1),enorm(1),vnorm(1))

      if (master) then
c        write(*,1001) rnorm(1),fnorm(1),enorm(1),vnorm(1)
 1001   format('*** Init. ||res||=',e10.4,' ||rhs||=',e10.4,
     &    ' ||err||=',e10.4,' ||var||=',e10.4)
      end if

C PROVENT THE CALCULATION OF rho FROM BREAK-DOWN
      if (rnorm(1) .le. small_number) then
        norm_last = 1.
      else
        norm_last = rnorm(1)
      end if
      nummg_l(1) = nummgits
      ilev       = 1

C-----------------------------------------------------------------------
C MULTIGRID SCHEDULE (ABLE TO GIVE V-, F- AND W-CYCLE)
C-----------------------------------------------------------------------

 10   if (nummg_l(1).gt.0 .and. enorm(1)/vnorm(1).gt.tolf) then

        goto_fine_grid = (ilev.eq.nummgl) .or. (nummg_l(ilev).eq.0)

C GOTO FINE GRID

        if (goto_fine_grid) then

C IF ON THE COARSEST GRID, SOLVE THE SYSTEM

          if (ilev .eq. nummgl) then
 
            call solversc(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,numcits,
     &           relax_pmgc,tolc,ra(jw1a),ia(iwk1a),rnorm(ilev),
     &           fnorm(ilev),enorm(ilev),vnorm(ilev))

            if (fcycle) mg_gamma = 1
 
          end if
 
C GOTO THE NEXT FINER GRID, DO PROLONGATION AND CORRECTION, AND THEN 
C DO POST-SMOOTHING

          if (ilev.gt.1 .and. nummgl.ne.1) then
 
            ilev = ilev - 1
            call prolongation(ilev,f,ra(jw1a),ia(iwk1a))

            call smooth(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,numpost,
     &           relax_pmg,tolc,ra(jw1a),ia(iwk1a),rnorm(ilev),
     &           fnorm(ilev),enorm(ilev),vnorm(ilev))

            if (ilev.eq.1 .and. fcycle) mg_gamma = 2
 
          end if
 
C COMPUTE THE REDUCTION FACTOR

          if (ilev .eq. 1) then
            rho_last  = rnorm(ilev)/norm_last
            norm_last = rnorm(ilev)
            if(master) then
c              write(*,1002) nummgits-nummg_l(ilev)+1,
c     &          rnorm(ilev),fnorm(ilev),enorm(ilev),vnorm(ilev),rho_last
 1002         format('*** i=',i3,' ||res||=',e10.4,' ||rhs||=',
     &          e10.4,' ||err||=',e10.4,' ||var||=',e10.4,' rho=',
     &          e10.4) 
            end if
          end if
 
          nummg_l(ilev) = nummg_l(ilev) - 1
 
C MOVE TO COARSER GRID

        else
 
C DO PRE-SMOOTHING

          call smooth(ilev,ap,aw,ae,as,an,al,ar,f,b,residu,numpre,
     &         relax_pmg,tolc,ra(jw1a),ia(iwk1a),rnorm(ilev),
     &         fnorm(ilev),enorm(ilev),vnorm(ilev))

          ilev = ilev + 1
          nummg_l(ilev) = mg_gamma
 
          call init_array(f,0.,1,ilev,ilev)
 
C COMPUTE THE RIGHT-HAND SIDE ON THE COARSE GRID

          call restriction(ilev,residu,b,ra(jw1a),ia(iwk1a))
 
        end if
 
C THIS goto STATEMENT CAN BE REPLACED EASILY BY do-while

        goto 10
 
      end if

!      if(master) then
!        write(*,1003)
!     &    nummgits-nummg_l(1),rnorm(1),fnorm(1),enorm(1),vnorm(1),
!     &    rho_last
! 1003  format('*** MG ',i4,
!     &   1x,'||res||=',e10.4,' ||rhs||=',e10.4, ' ||err||=',e10.4,
!     &   ' ||var||=',e10.4,' rho=',e10.4)
!      end if

 
      return
      end
