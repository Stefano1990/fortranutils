!-------------------------------
! Smagorinsky model subroutine
!-------------------------------

! Van Driest damping is done in virt_block

      subroutine subgrid_smagorinsky(gam,gamsc1,u,v,w,sc1,rho,
     &  vd_damping,dip,diu,djp,djv,dkp,dkw,
     &  hijp,hjip,hkip,hkjp,rp,bc,sc1var,sc1dis)

      implicit none

C PARAMETERS

#include "dim.inc"
#include "bc_type.inc"

C INCLUDE FILES

#include "cb.gridp.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "comm.inc"
#include "cb.comm.inc"

C Input/Output Parameters

      integer  bc(6)
      real     gam(nip1*njp1*nkp1), rho(nip1*njp1*nkp1)
      real     gamsc1(nip1*njp1*nkp1)
      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     sc1(nip1*njp1*nkp1)
      real     sc1var(nip1*njp1*nkp1),sc1dis(nip1*njp1*nkp1)
      real     vd_damping(nip1*njp1*nkp1)
      real     dip(nip1,njp1), diu(nip1,njp1), djp(nip1,njp1),
     &         djv(nip1,njp1)
      real     dkp(nkp1), dkw(nkp1)
      real     hijp(nip1,njp1), hjip(nip1,njp1), hkip(nip1,njp1),
     &         hkjp(nip1,njp1)
      real     rp(nip1,njp1)

C LOCAL

      integer  i, j, k, jk, ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr,
     &         ijkse, ijknw, ijkne, ijkel, ijknl, ijkwr, ijker, ijksr,
     &         ijknr
      integer  koff
      real     up, vp, wp
      real     di, dj, dk
      real     vwf, vef, vlf, vrf
      real     usf, unf, ulf, urf
      real     wwf, wef, wsf, wnf
      real     dudi, dudj, dudk, dvdi, dvdj, dvdk, dwdi, dwdj, dwdk
      real     radp, radw, rade, rads, radn
      real     wdrp, wdrw, wdre, wdrs, wdrn
      real     hki, hkj
      real     sijsij, delta_xyz, evisc
      real     dfdip, dfdjp, dfdkp, dfdijk
      real     fwf, fef, fsf, fnf, flf, frf
      real     sct, cssc1, cfv

C-----------------------------------------------------------------------
      data     sct/0.3/, cssc1/0.10/, cfv/0.1/
C-----------------------------------------------------------------------

C SOME CONSTANTS

      jo = nip1
      ko = nip1*njp1
      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C-----------------------------------------------------------------------
C LOOP OVER INTERIOR POINTS OF A BLOCK
C-----------------------------------------------------------------------

C THIS HAS NO WALL DAMPING!!!!!!!!!!!!!!
C - IF THE REYNOLDS NUMBER IS VERY HIGH, DAMPING IS NOT NECESSARY

      do k = 2,nk
        koff = ko*(k-1)
        do j = 2,nj
          jk = jo*(j-1) + koff
          do i = 2,ni

C????????????????????????????????????????????????????
C INITIALLY CODED FOR CLARITY - SPEED UP LATER
C????????????????????????????????????????????????????

            ijk   = i + jk

            ijkw  = ijk  - 1
            ijke  = ijk  + 1
            ijks  = ijk  - jo
            ijkn  = ijk  + jo
            ijkl  = ijk  - ko
            ijkr  = ijk  + ko

            ijkse = ijks + 1
            ijknw = ijkn - 1
            ijkne = ijkn + 1

            ijkel = ijke - ko
            ijknl = ijkn - ko
            ijkwr = ijkw + ko
            ijker = ijke + ko
            ijksr = ijks + ko
            ijknr = ijkn + ko

            fwf = 0.5*(sc1(ijk)+sc1(ijkw))
            fef = 0.5*(sc1(ijk)+sc1(ijke))
            fsf = 0.5*(sc1(ijk)+sc1(ijks))
            fnf = 0.5*(sc1(ijk)+sc1(ijkn))
            flf = 0.5*(sc1(ijk)+sc1(ijkl))
            frf = 0.5*(sc1(ijk)+sc1(ijkr))

            dfdip = (fef-fwf)/dip(i,j)
            dfdjp = (fnf-fsf)/djp(i,j)
            dfdkp = (frf-flf)/(rp(i,j)*dkw(k))

C CELL AVERAGE & CELL FACE VALUES            

            up = 0.5*(u(ijk) + u(ijke))
            vp = 0.5*(v(ijk) + v(ijkn))

            di = dip(i,j)
            dj = djp(i,j)

            vwf = 0.25*(v(ijk) + v(ijkn) + v(ijkw) + v(ijknw))
            vef = 0.25*(v(ijk) + v(ijkn) + v(ijke) + v(ijkne))

            usf = 0.25*(u(ijk) + u(ijke) + u(ijks) + u(ijkse))
            unf = 0.25*(u(ijk) + u(ijke) + u(ijkn) + u(ijkne))

            ulf = 0.25*(u(ijk) + u(ijke) + u(ijkl) + u(ijkel))
            urf = 0.25*(u(ijk) + u(ijke) + u(ijkr) + u(ijker))

            vlf = 0.25*(v(ijk) + v(ijkn) + v(ijkl) + v(ijknl))
            vrf = 0.25*(v(ijk) + v(ijkn) + v(ijkr) + v(ijknr))

C VELOCITY GRADIENTS : U & V IN I & J
        
            dudi = (u(ijke)-u(ijk))/di + vp*hijp(i,j)
            dvdi = (vef    -vwf   )/di - up*hijp(i,j)

            dudj = (unf    -usf   )/dj - vp*hjip(i,j)
            dvdj = (v(ijkn)-v(ijk))/dj + up*hjip(i,j)

C THE REST WITH RADIUS
C  - DIFFERENCE W/R => NO STRESS FOR SOLID BODY ROTATION ???????????????

           if (gotrad) then

              radp = rp(i  ,j  )
              radw = rp(i-1,j  )
              rade = rp(i+1,j  )
              rads = rp(i  ,j-1)
              radn = rp(i  ,j+1)

              dk = radp*dkp(k)

              hki  = hkip(i,j)
              hkj  = hkjp(i,j)

              wdrw = 0.5*(w(ijkw) + w(ijkwr))/radw
              wdre = 0.5*(w(ijke) + w(ijker))/rade
              wdrs = 0.5*(w(ijks) + w(ijksr))/rads
              wdrn = 0.5*(w(ijkn) + w(ijknr))/radn
              wdrp = 0.5*(w(ijk ) + w(ijkr ))/radp

              dwdi = radp*(wdre - wdrw)/(diu(i,j)+diu(i+1,j  ))
              dwdj = radp*(wdrn - wdrs)/(djv(i,j)+djv(i  ,j+1))

              wp = 0.5*(w(ijk) + w(ijkr))

C ?????? I THINK THE FOLLOWING TWO LINES ARE INCORRECT ???????
C             dudk = (urf    -ulf   )/dk
C             dvdk = (vrf    -vlf   )/dk
C ????????????????????????????????????????????????????????????
              dudk = (urf    -ulf   )/dk - wp*hki
              dvdk = (vrf    -vlf   )/dk - wp*hkj
              dwdk = (w(ijkr)-w(ijk))/dk + up*hki + vp*hkj

C THE REST WITHOUT RADIUS
C I THINK THIS IS NOT NECESSARY 

           else

             dk = dkp(k)

             wwf = 0.25*(w(ijk) + w(ijkr) + w(ijkw) + w(ijkwr))
             wef = 0.25*(w(ijk) + w(ijkr) + w(ijke) + w(ijker))
             wsf = 0.25*(w(ijk) + w(ijkr) + w(ijks) + w(ijksr))
             wnf = 0.25*(w(ijk) + w(ijkr) + w(ijkn) + w(ijknr))

             dwdi = (wef    -wwf   )/di
             dwdj = (wnf    -wsf   )/dj

             dudk = (urf    -ulf   )/dk
             dvdk = (vrf    -vlf   )/dk
             dwdk = (w(ijkr)-w(ijk))/dk

           end if

C ... SET GAM FOR OTHER DIRECTION & NOT OVERWRITE TW_SFK

            sijsij = dvdj*dvdj + dwdk*dwdk + dudi*dudi +
     &               0.5*((dudj+dvdi)*(dudj+dvdi) + 
     &                    (dvdk+dwdj)*(dvdk+dwdj) +
     &                    (dudk+dwdi)*(dudk+dwdi)   )

!            delta_xyz = (di*dj*dk)**.33333333
        delta_xyz = (di*dj*dk)**.33333333

           evisc = csles*delta_xyz*vd_damping(ijk)


            evisc = evisc*evisc*sqrt(2.0*sijsij)

            gam(ijk) = evisc
            if(csles .eq. 0.0) then
            gamsc1(ijk) = 0.0
            else 
            gamsc1(ijk) = gam(ijk)/sct
            endif  
C Temporary! Remove later

c            gam(ijk) = 0.0
c            gamsc1(ijk) = 0.0
            dfdijk = dfdip**2 + dfdjp**2 + dfdkp**2
            sc1var(ijk) = cfv*dfdijk*(delta_xyz**2)
            sc1dis(ijk) = sc1var(ijk) * sqrt(2.0*sijsij)

          end do
        end do
      end do

      return
      end
