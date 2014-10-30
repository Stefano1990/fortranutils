!-------------------------
! sigma subgrid scale model
!-------------------------

      subroutine subgrid_sigma(gam,gamsc1,u,v,w,sc1,rho,vd_damping,
     &  dip,diu,djp,djv,dkp,dkw,
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
      real     sct, cssc1, cfv,onethird
      real     g11,g22,g33,g12,g13,g21,g23,g31,g32,
     &         SdijSdij,Sdij,tracet
      real     invari1,invari2,invari3
      real     alpha1,alpha2,alpha3
      real     sigmafinal,sigma1,sigma2,sigma3
      real     pi,invari13rd,sqrtalpha
      real     gg11,gg22,gg33,tracegg

! What is sct?
      data     sct/0.3/, cssc1/0.10/, cfv/0.1/


! constants
      jo = nip1
      ko = nip1*njp1
      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1
      pi = 3.1415927
      onethird = 1.0/3.0 
          
      do k = 2,nk
        koff = ko*(k-1)
        do j = 2,nj
          jk = jo*(j-1) + koff
          do i = 2,ni

            ijk   = i + jk
            ! block sides W/E/S/...
            ijkw  = ijk  - 1
            ijke  = ijk  + 1
            ijks  = ijk  - jo
            ijkn  = ijk  + jo
            ijkl  = ijk  - ko
            ijkr  = ijk  + ko
            ! block edges (?)
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

!cell and cell face values

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

              dudk = (urf    -ulf   )/dk - wp*hki
              dvdk = (vrf    -vlf   )/dk - wp*hkj
              dwdk = (w(ijkr)-w(ijk))/dk + up*hki + vp*hkj

! if the domain is cartesian

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

!  gij**2 = gXX
!  gij**2 = gik*gkj
! redo these <<<<<<<<<<<<<<<<<<
! calculate capital G (name differnenty)
! use g11 et. with G to compute G**2
! G11 = g11g11(dudi etc.)+g21*g21+g31*g31 (double check)

            ! Do these remain the same?
! g is different for WALE model.
! These g calcs are for WALE model!
!            g11 = dudi*dudi + dudj*dvdi + dudk*dwdi
!            g22 = dvdi*dudj + dvdj*dvdj + dvdk*dwdj
!            g33 = dwdi*dudk + dwdj*dvdk + dwdk*dwdk
!            g12 = dudi*dudj + dudj*dvdj + dudk*dwdj
!            g13 = dudi*dudk + dudj*dvdk + dudk*dwdk
!            g21 = dvdi*dudi + dvdj*dvdi + dvdk*dwdi 
!            g23 = dvdi*dudk + dvdj*dvdk + dvdk*dwdk
!            g31 = dwdi*dudi + dwdj*dvdi + dwdk*dwdi
!            g32 = dwdi*dudj + dwdj*dvdj + dwdk*dwdj
! End of WALE g

! Sigma gs
        g11 = dudi*dudi + dvdi*dvdi + dwdi*dwdi
        g12 = dudi*dudj + dvdi*dvdj + dwdi*dwdj
        g13 = dudi*dudk + dvdi*dvdk + dwdi*dwdk
        g21 = dudj*dudi + dvdj*dvdi + dwdj*dwdi
        g22 = dudj*dudj + dvdj*dvdj + dwdj*dwdj
        g23 = dudj*dudk + dvdj*dvdk + dwdj*dwdk
        g31 = dudk*dudi + dvdk*dvdi + dwdk*dwdi
        g32 = dudk*dudj + dvdk*dvdj + dwdk*dwdj
        g33 = dudk*dudk + dvdk*dvdk + dwdk*dwdk

! Calculate tr(G**2)
        gg11 = g11*g11 + g12*g21 + g13*g31
        gg22 = g21*g12 + g22*g22 + g23*g32
        gg33 = g31*g13 + g32*g23 + g33*g33
        tracegg = gg11 + gg22 + gg33
! 2. compute invariants 
! 2.1 Trace of G
        invari1 = g11 + g22 + g33
! 2.2 Invariant 2
! Fix this.
        invari2 = 0.5*(invari1**2-tracegg)
! 2.3 Determinant of G
        invari3 =  g11*(g22*g33 - g23*g32)
     &            - g12*(g21*g33 - g23*g31)
     &            + g13*(g21*g32 - g22*g31)
! 3. Compute angles from invariants
        alpha1 = (invari1**2*0.11111111) - (invari2*onethird)
        alpha2 = (invari1**3*0.03703703703) 
     &             - (invari1*invari2*0.16666666666)
     &             + (invari3*0.5)
        alpha3 = onethird*acos(alpha2/(alpha1**1.5))

! 4. Compute singular values
! Should I cache sqrt(alpha1) so I compute it 1 time instead of 3?
! Check speed difference for this.
! Can also cache invari1/3.0
        invari13rd = invari1*onethird
        sqrtalpha = sqrt(alpha)
        sigma1 = (invari13rd+2*sqrtalpha*cos(alpha3))**0.5
        sigma2 = (invari13rd-2*sqrtalpha*cos(pi*onethird+alpha3))**0.5
        sigma3 = (invari13rd-2*sqrtalpha*cos(pi*onethird-alpha3))**0.5

! 5. Compute sigma.
        sigmafinal = (sigma3*(sigma1-sigma2)*(sigma2-sigma3))
     &                  /(sigma1**2.0)
        if(isnan(sigmafinal)) then
!         write(*,*) 'invari1:',invari1,'invari2:',invari2,'in3:',invari3
!         write(*,*) 'inside:',alpha2/(alpha1**1.5)
!         write(*,*) 'acos:',acos(alpha2/(alpha1**1.5))
        if(pe .eq. 6) then 
        write(*,*) i,j,k,ijk,pe
        endif  
         !write(*,*) '...'
!         STOP
        endif
! not needed.
! trace part - dij gkk

!           tracet = Onethird*(g11+g22+g33)            

! not needed.
!Sdij = 0.5*(gij**2 + gji**2) - (1/3)dij gkk
!         Sdij = (g11 - tracet) + (g22 - tracet) + (g33 - tracet)
!     &        + 0.5*(g12+g21) + 0.5*(g13+g31) + 0.5*(g23+g32)
!         SdijSdij = Sdij**2 
!
        delta_xyz = (di*dj*dk)**.33333333
!        write(*,*) 'delta_xyz:',delta_xyz

!Cs here refers to the wale constant
        evisc = (csles*delta_xyz)**2
!        if(isnan(evisc)) write(*,*) 'evisc:',evisc

!        gam(ijk) = evisc*( (SdijSdij)**1.5 / 
!     &          (sijsij**2.5 + SdijSdij**1.25))
        gam(ijk) = evisc*sigmafinal   
!        write(*,*) 'gam(ijk):',gam(ijk)
            if(csles .eq. 0.0) then
            gamsc1(ijk) = 0.0
            else 
            gamsc1(ijk) = gam(ijk)/sct
            endif  

            dfdijk = dfdip**2 + dfdjp**2 + dfdkp**2
            sc1var(ijk) = cfv*dfdijk*(delta_xyz**2)
            sc1dis(ijk) = sc1var(ijk) * sqrt(2.0*sijsij)

          end do
        end do
      end do

      return
      end
