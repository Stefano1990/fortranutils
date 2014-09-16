C=======================================================================
C SET BOUNDARY VALUES FOR VIRTUAL BOUNDARY POINTS
C=======================================================================

C NOTE:
C 1. NOT FOR GENERAL BOUNDARY CONDITIONS IN WHICH value IS A FUNCTION
C    OF GRID LEVEL, FACES AND LOCATIONS
C 2. ARRAY gam IS USED TO PASSED THE VAN DRIEST DAMPING FACTOR TO
C    SUBROUTNE set_gam_block


      subroutine set_virt_bc_block(iva,nip1,njp1,nkp1,
     &  f,fm,rho,gam,vd_damping,value,
     &  iside,k,k1,io,jo,ko,kko,
     &  kappa,visc,aplus,ylaminar,
     &  walls,init_vdd,twf,
     &  dip,diu,div,dix,djp,dju,djv,djx,ru,rv,dkp,dkw)

      implicit none

#include "les.inc"

      integer iva, iside, k, k1, kk
      integer nip1, njp1, nkp1
      integer io, jo, ko, kko, kko2
      logical walls(6), init_vdd(6)

      real    f(nip1*njp1*nkp1), fm(nip1*njp1*nkp1)
      real    rho(nip1*njp1*nkp1), gam(nip1*njp1*nkp1)
      real    vd_damping(nip1*njp1*nkp1)
      real    value, kappa, visc, ylaminar, aplus
      real    twf(*)
      real    dip(*), diu(*), div(*), dix(*),
     &        djp(*), dju(*), djv(*), djx(*),
     &        ru(*) , rv(*) , dkp(*), dkw(*)

C LOCAL

      integer i, j, ij, ijk, ijk_next
      integer jkoff, koff, joff1, koff1, koff0, koff2
      integer kvd, kref, kstart, kfinish, kstep, kstr

      logical wall, wallwf, wallef, wallsf, wallnf, walllf, wallrf

      real    dist, tw_mean, ut_mean, yplus_mean, tw
      real    wewall, snwall, lrwall
      real    rhof, rhoflast, rhofnlast


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      wallwf = walls(1)
      wallef = walls(2)
      wallsf = walls(3)
      wallnf = walls(4)
      walllf = walls(5)
      wallrf = walls(6)
      wall = wallwf .or. wallef .or.
     &       wallsf .or. wallnf .or.
     &       walllf .or. wallrf
      wewall = 1.
      snwall = 1.
      lrwall = 1.
      if (wallef) wewall = -1.
      if (wallnf) snwall = -1.
      if (wallrf) lrwall = -1.

      koff = 1 + ko*(k-1)
      kko2 = kko + kko

C IF WALL AND VELOCITY, SET THE VIRTUAL VELOCITY ACCORDING TO THE
C LOG-LAW WHEN yplus EXCEEDS ylaminar

      if (wall .and. iva.ge.SV_U .and. iva.le.SV_W) then
      
C-----------------------------------------------------------------------
C u-VELOCITY
C-----------------------------------------------------------------------

       if (iva .eq. SV_U) then

C - SOUTH WALL OR NORTH WALL

          if (wallsf .or. wallnf) then

            if (k1 .eq. 2) then
               kstart = k1
               kfinish = k1 + nkp1/2
c               kfinish = k1
c               kfinish = k1 + 25
               kstep = 1
               kstr = k1 + 1
            else
               kstart = k1
               kfinish = k1 - nkp1/2
c               kfinish = k1
c               kfinish = k1 - 25
               kstep = -1
               kstr = k1 - 1
            endif
            koff0 = nip1*(k1-1)

            if (init_vdd(iside)) then
              init_vdd(iside) = .false.
              do kvd = kstart, kfinish, kstep

                 if (k1 .eq. 2) then
                    kref = kvd - 1
                 else
                    kref = kvd
                 endif

                 koff2 = 1 + ko*(kref - 1)
                 do j = 1,njp1
                   jkoff = jo*(j-1) + koff2
                   do i = 1,nip1
                      ijk = io*(i-1) + jkoff
                      ijk_next = ijk + kko
                      vd_damping(ijk_next) = 0.
                   end do
                 end do
              end do
            end if

            do kvd = kstart, kfinish, kstep

            if (k1 .eq. 2) then
               kref = kvd - 1
            else
               kref = kvd
            endif

            koff1 = nip1*(kvd-1)
            koff2 = 1 + ko*(kref - 1)
            do j = 1,njp1
              jkoff = jo*(j-1) + koff2
              rhoflast = rho(jkoff)
              rhofnlast = rho(jkoff + kko)
              do i = 1,nip1
                ijk      = io*(i-1) + jkoff
                ijk_next = ijk + kko

                rhof       = 0.25*(rho(ijk   )+rho(ijk_next   )
     &                           +rho(ijk-io)+rho(ijk_next-io))
c     &                            +rhoflast   +rhofnlast)
                rhoflast = rho(ijk)
                rhofnlast = rho(ijk_next)
                tw         = twf(j+njp1*(i-1))
                tw_mean    = abs(tw)
                ut_mean    = sqrt(tw_mean/rhof)
                if (kvd .eq. kstart) then
                   dist       = 0.5*dju(i+koff0)
                else
                   dist       = 0.5*dju(i+koff0)
                   do kk = kstr, kvd, kstep
                      dist = dist + djx(i+nip1*(kk-1))
                   enddo
c                   dist = 0.5*dju(i+koff0)+abs(ru(i+koff1)-ru(i+koff0))
                endif
                yplus_mean = dist*ut_mean*rhof/visc

                vd_damping(ijk_next) = 
c     &            max(1. - exp(-yplus_mean/aplus),vd_damping(ijk_next))
c     &            max(1. - exp(-yplus_mean/aplus),0.0)
     &            amin1(1. - exp(-yplus_mean/aplus), 1.0)
                if (kvd .eq. kstart) then
                   if (yplus_mean.gt.ylaminar) then
                     f(ijk) = f(ijk_next)
     &                 - snwall*sign(1.,tw)*ut_mean*djx(i+koff0)
     &                   /(kappa*dist)
                   else
                     f(ijk) = 2.*value - f(ijk_next)
                   end if
                endif

              end do
            end do
            end do

          end if

C - LEFT WALL OR RIGHT WALL

          if (walllf .or. wallrf) then

            if (init_vdd(iside)) then
              init_vdd(iside) = .false.
              do j = 1,njp1
                jkoff = jo*(j-1) + koff + kko
                do i = 1,nip1
                  vd_damping(io*(i-1) + jkoff) = 0.
                end do
              end do
            end if
            do j = 1,njp1
              jkoff = jo*(j-1) + koff
              joff1 = nip1*(j-1)
              rhoflast = rho(jkoff)
              rhofnlast = rho(jkoff + kko)
              do i = 1,nip1
                ijk      = io*(i-1) + jkoff
                ijk_next = ijk + kko
                ij       = i + joff1

                rhof       = 0.25*(rho(ijk   )+rho(ijk_next   )
c    &                            +rho(ijk-io)+rho(ijk_next-io))
     &                            +rhoflast   +rhofnlast)
                rhoflast = rho(ijk)
                rhofnlast = rho(ijk_next)
                dist       = 0.5*ru(ij)*dkp(k1)
                tw         = twf(ij)
                tw_mean    = abs(tw)
                ut_mean    = sqrt(tw_mean/rhof)
                yplus_mean = dist*ut_mean*rhof/visc

                vd_damping(ijk_next) = 
     &            max(1. - exp(-yplus_mean/aplus),vd_damping(ijk_next))
                if (yplus_mean.gt.ylaminar) then
                  f(ijk) = f(ijk_next)
     &              - lrwall*sign(1.,tw)*ut_mean*ru(ij)*dkw(k1)
     &                /(kappa*dist)
                else
                  f(ijk) = 2.*value - f(ijk_next)
                end if

              end do
            end do

          end if

C-----------------------------------------------------------------------
C v-VELOCITY
C-----------------------------------------------------------------------
        else if (iva .eq. SV_V) then

C - WEST WALL OR EAST WALL

          if (wallwf .or. wallef) then

            if (init_vdd(iside)) then
              init_vdd(iside) = .false.
              do j = 1,njp1
                jkoff = jo*(j-1) + koff + kko
                do i = 1,nip1
                  vd_damping(io*(i-1) + jkoff) = 0.
                end do
              end do
            end if
            do j = 1,njp1
              jkoff = jo*(j-1) + koff
              rhoflast = rho(jkoff)
              rhofnlast = rho(jkoff + kko)
              do i = 1,nip1
                ijk      = io*(i-1) + jkoff
                ijk_next = ijk + kko

                rhof       = 0.25*(rho(ijk   )+rho(ijk_next   )
c    &                            +rho(ijk-io)+rho(ijk_next-io))
     &                            +rhoflast   +rhofnlast)
                rhoflast = rho(ijk)
                rhofnlast = rho(ijk_next)
                dist       = 0.5*div(k1+nkp1*(i-1))
                tw         = twf(j+njp1*(i-1))
                tw_mean    = abs(tw)
                ut_mean    = sqrt(tw_mean/rhof)
                yplus_mean = dist*ut_mean*rhof/visc

                vd_damping(ijk_next) = 
     &            max(1. - exp(-yplus_mean/aplus),vd_damping(ijk_next))
                if (yplus_mean.gt.ylaminar) then
                  f(ijk) = f(ijk_next)
     &             - wewall*sign(1.,tw)*ut_mean*dix(k1+nkp1*(i-1))
     &               /(kappa*dist)
                else
                  f(ijk) = 2.*value - f(ijk_next)
                end if

              end do
            end do

          end if

C - LEFT WALL OR RIGHT WALL

          if (walllf .or. wallrf) then

            if (init_vdd(iside)) then
              init_vdd(iside) = .false.
              do j = 1,njp1
                jkoff = jo*(j-1) + koff + kko
                do i = 1,nip1
                  vd_damping(io*(i-1) + jkoff) = 0.
                end do
              end do
            end if
            do i = 1,nip1
              rhoflast = rho(koff)
              rhofnlast = rho(koff + kko)
              do j = 1,njp1
                jkoff = jo*(j-1) + koff
                joff1 = nip1*(j-1)

                ijk      = io*(i-1) + jkoff
                ijk_next = ijk + kko
                ij       = i + joff1

                rhof       = 0.25*(rho(ijk   )+rho(ijk_next   )
c    &                            +rho(ijk-jo)+rho(ijk_next-jo))
     &                            +rhoflast   +rhofnlast)
                rhoflast = rho(ijk)
                rhofnlast = rho(ijk_next)
                dist       = 0.5*rv(ij)*dkp(k1)
                tw         = twf(ij)
                tw_mean    = abs(tw)
                ut_mean    = sqrt(tw_mean/rhof)
                yplus_mean = dist*ut_mean*rhof/visc

                vd_damping(ijk_next) = 
     &            max(1. - exp(-yplus_mean/aplus),vd_damping(ijk_next))
                if (yplus_mean.gt.ylaminar) then
                  f(ijk) = f(ijk_next)
     &              - lrwall*sign(1.,tw)*ut_mean*rv(ij)*dkw(k1)
     &                /(kappa*dist)
                else
                  f(ijk) = 2.*value - f(ijk_next)
                end if

              end do
            end do

          end if

C-----------------------------------------------------------------------
C w-VELOCITY
C-----------------------------------------------------------------------
        else

C WEST WALL OR EAST WALL

          if (wallwf .or. wallef) then

            if (init_vdd(iside)) then
              init_vdd(iside) = .false.
              do j = 1,njp1
                jkoff = jo*(j-1) + koff + kko
                do i = 1,nip1
                  vd_damping(io*(i-1) + jkoff) = 0.
                end do
              end do
            end if
            do i = 1,nip1
              rhoflast = rho(koff)
              rhofnlast = rho(koff + kko)
              do j = 1,njp1
                jkoff = jo*(j-1) + koff

                ijk      = io*(i-1) + jkoff
                ijk_next = ijk + kko

                rhof       = 0.25*(rho(ijk   )+rho(ijk_next   )
c    &                            +rho(ijk-jo)+rho(ijk_next-jo))
     &                            +rhoflast   +rhofnlast)
                rhoflast = rho(ijk)
                rhofnlast = rho(ijk_next)
                dist       = 0.5*dip(k1+nkp1*(i-1))
                tw         = twf(j+njp1*(i-1))
                tw_mean    = abs(tw)
                ut_mean    = sqrt(tw_mean/rhof)
                yplus_mean = dist*ut_mean*rhof/visc

                vd_damping(ijk_next) = 
     &            max(1. - exp(-yplus_mean/aplus),vd_damping(ijk_next))
                if (yplus_mean.gt.ylaminar) then
                  f(ijk) = f(ijk_next)
     &             - wewall*sign(1.,tw)*ut_mean*diu(k1+nkp1*(i-1))
     &               /(kappa*dist)
                else
                  f(ijk) = 2.*value - f(ijk_next)
                end if

              end do
            end do

          end if

C SOUTH WALL OR NORTH WALL

          if (wallsf .or. wallnf) then

            if (init_vdd(iside)) then
              init_vdd(iside) = .false.
              do j = 1,njp1
                jkoff = jo*(j-1) + koff + kko
                do i = 1,nip1
                  vd_damping(io*(i-1) + jkoff) = 0.
                end do
              end do
            end if
            koff1 = nip1*(k1-1)
            do i = 1,nip1
              rhoflast = rho(koff)
              rhofnlast = rho(koff + kko)
              do j = 1,njp1
                jkoff = jo*(j-1) + koff

                ijk      = io*(i-1) + jkoff
                ijk_next = ijk + kko

                rhof       = 0.25*(rho(ijk   )+rho(ijk_next   )
c    &                            +rho(ijk-jo)+rho(ijk_next-jo))
     &                            +rhoflast   +rhofnlast)
                rhoflast = rho(ijk)
                rhofnlast = rho(ijk_next)
                dist       = 0.5*djp(i+koff1)
                tw         = twf(j+njp1*(i-1))
                tw_mean    = abs(tw)
                ut_mean    = sqrt(tw_mean/rhof)
                yplus_mean = dist*ut_mean*rhof/visc

                vd_damping(ijk_next) = 
     &            max(1. - exp(-yplus_mean/aplus),vd_damping(ijk_next))
                if (yplus_mean.gt.ylaminar) then
                  f(ijk) = f(ijk_next)
     &              - snwall*sign(1.,tw)*ut_mean*djv(i+koff1)/
     &                (kappa*dist)
                else
                  f(ijk) = 2.*value - f(ijk_next)
                end if

              end do
            end do

          end if

        end if

C-----------------------------------------------------------------------
C OTHER VARIALES, OR VELOCITIES ON NON-WALL BOUNDARY
C-----------------------------------------------------------------------
      else
       write(*,*) 'in here doing stuff'
       
        koff = 1 + ko*(k-1)
        do j = 1,njp1
          jkoff = jo  *(j-1) + koff
          do i = 1,nip1
            ijk = io*(i-1) + jkoff
            f(ijk) = 2.*value - f(ijk+kko)
          end do
        end do

      end if

      return
      end
