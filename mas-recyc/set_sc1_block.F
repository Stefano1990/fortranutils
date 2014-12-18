C=======================================================================
C SET SCALAR FIELD
C=======================================================================

      subroutine set_sc1_block(nip1,njp1,nkp1,io,jo,ko,dtime,
     &                     dip, diu, rdiv,
     &                     djp, djv, rdju,
     &                     didjp, rdidjp,rp,
     &                     dkp, dkw,
     &                     u,v,w,p,sc1,hsc1,rho,gamsc1,gamlam,visc,
     &                     kstart,iststep,iblk)

      implicit none

#include "dim.inc"
#include "const.inc"
#include "bc_type.inc"
#include "cb.lgrid.inc"
#include "tvdgrad.inc"

C PARAMETERS

C ARGUMENTS

      integer      nip1, njp1, nkp1
      integer      io, jo, ko
      integer      kstart
      integer      iblk
      real         dtime

      real         dip(nip1,njp1), diu(nip1,njp1), rdiv(nip1,njp1)
      real         djp(nip1,njp1), djv(nip1,njp1), rdju(nip1,njp1)
      real         didjp(nip1,njp1), rdidjp(nip1,njp1)
      real         rp(nip1,njp1)
      real         dkp(nkp1), dkw(nkp1)

      real         u(nip1*njp1*nkp1)
      real         v(nip1*njp1*nkp1)
      real         w(nip1*njp1*nkp1)
      real         p(nip1*njp1*nkp1)
      real         sc1(nip1*njp1*nkp1)
      real         hsc1(nip1*njp1*nkp1)
      real         rho(nip1*njp1*nkp1)
      real         gamsc1(nip1*njp1*nkp1),gamlam(nip1*njp1*nkp1)
      real         visc

      integer      iststep

C LOCAL

      integer      i, j, k
      integer      ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr, ijk2
      integer      ijkww, ijkss, ijkll
      integer      ijkee, ijknn, ijkrr
      integer      ijksw, ijknw, ijkwl, ijkwr
      integer      iblk_next, iside_next
      integer      ij, ik, jk
      integer      ni, nj, nk
      integer      bc
      
      real         rhowf, rhoef, rhosf, rhonf, rholf, rhorf, rhop
      real         gamwf, gamef, gamsf, gamnf, gamlf, gamrf
      real         uwf, uef, vsf, vnf, wlf, wrf
      real         dfdxwf, dfdxef, dfdxwwf,dfdxeef
      real         dfdysf, dfdynf, dfdyssf,dfdynnf
      real         dfdzlf, dfdzrf, dfdzllf,dfdzrrf
      real         rwst, rest, rsth, rnth, rlft, rrht
      real         limw, lime, lims, limn, liml, limr
      real         fwf, fef, fsf, fnf, flf, frf
      real         areawf, areaef, areasf, areanf, areak
      real         gwf, gef, gsf, gnf, glf, grf
      real         galwf, galef, galsf, galnf, gallf, galrf
      real         twf, tef, tsf, tnf, tlf, trf
      real         h
      real         st1, st2, stden
      real         vol
      real         convf, difff
      real         dwdy, dvdz, dwdx, dudz
      real         lambda, sclam
      real         sc1tmp

C DATA

      data lambda/4.0/, sclam/0.7/

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------

      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

      
C-----------------------------------------------------------------------
C SET SC1 & HSC1
C-----------------------------------------------------------------------

      do j = 2, nj
        do i = 2, ni

C SET I,J METRICS

          do k = 2, nk

            ijk  = i + jo*(j-1) + ko*(k-1)

            ijkw = ijk - io
            ijke = ijk + io
            ijks = ijk - jo
            ijkn = ijk + jo
            ijkl = ijk - ko
            ijkr = ijk + ko

            if (i .gt. 2) then
               ijkww = ijkw - io
            endif

            if (i .lt. ni) then
               ijkee = ijke + io
            endif

            if (j .gt. 2) then
               ijkss = ijks - jo
            endif

            if (j .lt. nj) then
               ijknn = ijkn + jo
            endif

            if (k .gt. 2) then
               ijkll = ijkl - ko
            endif

            if (k .lt. nk) then
               ijkrr = ijkr + ko
            endif

C-----------------------------------------------------------------------
C FACE AND CELL AVERAGE VALUES
C-----------------------------------------------------------------------

            rhowf = 0.5*(rho(ijk)+rho(ijkw))
            rhoef = 0.5*(rho(ijk)+rho(ijke))
            rhosf = 0.5*(rho(ijk)+rho(ijks))
            rhonf = 0.5*(rho(ijk)+rho(ijkn))
            rholf = 0.5*(rho(ijk)+rho(ijkl))
            rhorf = 0.5*(rho(ijk)+rho(ijkr))
            rhop = rho(ijk)

            gamwf = 0.5*(gamsc1(ijk)+gamsc1(ijkw))*rhowf + 
     &              0.5*(gamlam(ijk)+gamlam(ijkw))
            gamef = 0.5*(gamsc1(ijk)+gamsc1(ijke))*rhoef + 
     &              0.5*(gamlam(ijk)+gamlam(ijke))
            gamsf = 0.5*(gamsc1(ijk)+gamsc1(ijks))*rhosf + 
     &              0.5*(gamlam(ijk)+gamlam(ijks))
            gamnf = 0.5*(gamsc1(ijk)+gamsc1(ijkn))*rhonf + 
     &              0.5*(gamlam(ijk)+gamlam(ijkn))
            gamlf = 0.5*(gamsc1(ijk)+gamsc1(ijkl))*rholf + 
     &              0.5*(gamlam(ijk)+gamlam(ijkl))
            gamrf = 0.5*(gamsc1(ijk)+gamsc1(ijkr))*rhorf + 
     &              0.5*(gamlam(ijk)+gamlam(ijkr))

            uwf = u(ijk)
            uef = u(ijke)
            vsf = v(ijk)
            vnf = v(ijkn)
            wlf = w(ijk)
            wrf = w(ijkr)

c            if (j .eq. 2) write(192,*) "j=",j," vsf=",vsf," vnf=",vnf
c            if (j .eq. nj) write(192,*) "j=",j," vsf=",vsf," vnf=",vnf

            dfdxwf = (sc1(ijk)-sc1(ijkw))/diu(i,j)
            dfdxef = (sc1(ijke)-sc1(ijk))/diu(i+1,j)
            if (i .eq. 2) then
               bc = bl_bc(1,iblk)
               if (bc .ne. BCT_INTERFACE) then
                  dfdxwwf = dfdxwf
               else
                  iblk_next = bl2bl(1,iblk)
                  iside_next = bl2side(1,iblk)
                  dfdxwwf = dfdxe(iblk_next,j,k)
               endif
            else
               dfdxwwf = (sc1(ijkw)-sc1(ijkww))/diu(i-1,j)
            endif
            if (i .eq. ni) then
               bc = bl_bc(2,iblk)
               if (bc .ne. BCT_INTERFACE) then
                  dfdxeef = dfdxef
               else
                  iblk_next = bl2bl(2,iblk)
                  iside_next = bl2side(2,iblk)
                  dfdxeef = dfdxw(iblk_next,j,k)
               endif
            else
               dfdxeef = (sc1(ijkee)-sc1(ijke))/diu(i+2,j)
            endif

            dfdysf = (sc1(ijk)-sc1(ijks))/djv(i,j)
            dfdynf = (sc1(ijkn)-sc1(ijk))/djv(i,j+1)
            if (j .eq. 2) then
               bc = bl_bc(3,iblk)
               if (bc .ne. BCT_INTERFACE) then
                  dfdyssf = dfdysf
               else
                  iblk_next = bl2bl(3,iblk)
                  iside_next = bl2side(3,iblk)
                  dfdyssf = dfdyn(iblk_next,i,k)
               endif
            else
               dfdyssf = (sc1(ijks)-sc1(ijkss))/djv(i,j-1)
            endif
            if (j .eq. nj) then
               bc = bl_bc(4,iblk)
               if (bc .ne. BCT_INTERFACE) then
                  dfdynnf = dfdynf
               else
                  iblk_next = bl2bl(4,iblk)
                  iside_next = bl2side(4,iblk)
                  dfdynnf = dfdys(iblk_next,i,k)
               endif
            else
               dfdynnf = (sc1(ijknn)-sc1(ijkn))/djv(i,j+2)
            endif

c            if (i .eq. 36 .and. k .eq. 29) then
c               write (192,*) j, djv(i,j),djv(i,j-1),djv(i,j+1)
c            endif
c            dfdzlf = (sc1(ijk)-sc1(ijkl))/(rp(i,j)*dkp(k))
c            dfdzrf = (sc1(ijkr)-sc1(ijk))/(rp(i,j)*dkp(k+1))
            dfdzlf = (sc1(ijk)-sc1(ijkl))/(rp(i,j)*dkp(k-1))
            dfdzrf = (sc1(ijkr)-sc1(ijk))/(rp(i,j)*dkp(k))
            if (k .eq. 2) then
               bc = bl_bc(5,iblk)
               if (bc .ne. BCT_INTERFACE) then
                  dfdzllf = dfdzlf
               else
                  iblk_next = bl2bl(5,iblk)
                  iside_next = bl2side(5,iblk)
                  dfdzllf = dfdzr(iblk_next,i,j)
               endif
            else
c               dfdzllf = (sc1(ijkl)-sc1(ijkll))/(rp(i,j)*dkp(k-1))
               dfdzllf = (sc1(ijkl)-sc1(ijkll))/(rp(i,j)*dkp(k-2))
            endif
            if (k .eq. nk) then
               bc = bl_bc(6,iblk)
               if (bc .ne. BCT_INTERFACE) then
                  dfdzrrf = dfdzrf
               else
                  iblk_next = bl2bl(6,iblk)
                  iside_next = bl2side(6,iblk)
                  dfdzrrf = dfdzl(iblk_next,i,j)
               endif
            else
c               dfdzrrf = (sc1(ijkrr)-sc1(ijkr))/(rp(i,j)*dkp(k+2))
               dfdzrrf = (sc1(ijkrr)-sc1(ijkr))/(rp(i,j)*dkp(k+1))
            endif

C CONSECUTIVE GRADIENT PARAMETER
 
            if (uwf .ge. 0.0) then
               rwst = (dfdxwwf + small_number)/(dfdxwf + small_number)
            else
               rwst = (dfdxef + small_number)/(dfdxwf + small_number)
            endif
            if (uef .gt. 0.0) then
               rest = (dfdxwf + small_number)/(dfdxef + small_number)
            else
               rest = (dfdxeef + small_number)/(dfdxef + small_number)
            endif
            if (vsf .ge. 0.0) then
               rsth = (dfdyssf + small_number)/(dfdysf + small_number)
            else
               rsth = (dfdynf + small_number)/(dfdysf + small_number)
            endif
            if (vnf .gt. 0.0) then
               rnth = (dfdysf + small_number)/(dfdynf + small_number)
            else
               rnth = (dfdynnf + small_number)/(dfdynf + small_number)
            endif
            if (wlf .ge. 0.0) then
               rlft = (dfdzllf + small_number)/(dfdzlf + small_number)
            else
               rlft = (dfdzrf + small_number)/(dfdzlf + small_number)
            endif
            if (wrf .gt. 0.0) then
               rrht = (dfdzlf + small_number)/(dfdzrf + small_number)
            else
               rrht = (dfdzrrf + small_number)/(dfdzrf + small_number)
            endif

C TVD FLUX LIMITERS: ROE'S ORIGINAL, "minmod"

            limw = amax1(0.0,amin1(lambda*rwst,1.0))
            lime = amax1(0.0,amin1(lambda*rest,1.0))
            lims = amax1(0.0,amin1(lambda*rsth,1.0))
            limn = amax1(0.0,amin1(lambda*rnth,1.0))
            liml = amax1(0.0,amin1(lambda*rlft,1.0))
            limr = amax1(0.0,amin1(lambda*rrht,1.0))

C Remove comments for upwind differencing

c            limw = 0.0
c            lime = 0.0
c            lims = 0.0
c            limn = 0.0
c            liml = 0.0
c            limr = 0.0

C Remove comments for central differencing

c            limw = 1.0
c            lime = 1.0
c            lims = 1.0
c            limn = 1.0
c            liml = 1.0
c            limr = 1.0

            if (uwf .gt. 0.0) then
               fwf = sc1(ijkw) + 0.5*limw*(sc1(ijk)-sc1(ijkw))
            else if (uwf .lt. 0.0) then
               fwf = sc1(ijk) + 0.5*limw*(sc1(ijkw)-sc1(ijk))
            else
               fwf = 0.5*(sc1(ijk) + sc1(ijkw))
            endif
            if (uef .gt. 0.0) then
               fef = sc1(ijk) + 0.5*lime*(sc1(ijke)-sc1(ijk))
            else if (uef .lt. 0.0) then
               fef = sc1(ijke) + 0.5*lime*(sc1(ijk)-sc1(ijke))
            else
               fef = 0.5*(sc1(ijk) + sc1(ijke))
            endif
            if (vsf .gt. 0.0) then
               fsf = sc1(ijks) + 0.5*lims*(sc1(ijk)-sc1(ijks))
            else if (vsf .lt. 0.0) then
               fsf = sc1(ijk) + 0.5*lims*(sc1(ijks)-sc1(ijk))
            else
               fsf = 0.5*(sc1(ijk) + sc1(ijks))
            endif
            if (vnf .gt. 0.0) then
               fnf = sc1(ijk) + 0.5*limn*(sc1(ijkn)-sc1(ijk))
            else if (vnf .lt. 0.0) then
               fnf = sc1(ijkn) + 0.5*limn*(sc1(ijk)-sc1(ijkn))
            else
               fnf = 0.5*(sc1(ijk) + sc1(ijkn))
            endif
            if (wlf .gt. 0.0) then
               flf = sc1(ijkl) + 0.5*liml*(sc1(ijk)-sc1(ijkl))
            else if (wlf .lt. 0.0) then
               flf = sc1(ijk) + 0.5*liml*(sc1(ijkl)-sc1(ijk))
            else
               flf = 0.5*(sc1(ijk) + sc1(ijkl))
            endif
            if (wrf .gt. 0.0) then
               frf = sc1(ijk) + 0.5*limr*(sc1(ijkr)-sc1(ijk))
            else if (wrf .lt. 0.0) then
               frf = sc1(ijkr) + 0.5*limr*(sc1(ijk)-sc1(ijkr))
            else
               frf = 0.5*(sc1(ijk) + sc1(ijkr))
            endif

            if (i .eq. 2) then
c               fwf = 0.5*(sc1(ijk) + sc1(ijkw))
c               fef = 0.5*(sc1(ijk) + sc1(ijke))
            endif
            if (i .eq. ni) then
c               fwf = 0.5*(sc1(ijk) + sc1(ijkw))
c               fef = 0.5*(sc1(ijk) + sc1(ijke))
            endif
            if (j .eq. 2) then
c               fsf = 0.5*(sc1(ijk) + sc1(ijks))
c               fnf = 0.5*(sc1(ijk) + sc1(ijkn))
            endif
            if (j .eq. nj) then
c               fsf = 0.5*(sc1(ijk) + sc1(ijks))
c               fnf = 0.5*(sc1(ijk) + sc1(ijkn))
            endif
            if (k .eq. 2) then
c               flf = 0.5*(sc1(ijk) + sc1(ijkl))
c               frf = 0.5*(sc1(ijk) + sc1(ijkr))
            endif
            if (k .eq. nk) then
c               flf = 0.5*(sc1(ijk) + sc1(ijkl))
c               frf = 0.5*(sc1(ijk) + sc1(ijkr))
            endif

c            if (i .eq. 2) then
c               if (uwf .ge. 0.0) then
c                  fwf = 0.5*(sc1(ijk) + sc1(ijkw))
c               else
c                  fwf = sc1(ijk)
c               endif
c            endif
c
c            if (i .eq. ni) then
c               if (uef .ge. 0.0) then
c                  fef = sc1(ijk)
c               else
c                  fef = 0.5*(sc1(ijk) + sc1(ijke))
c               endif
c            endif
c
c            if (j .eq. 2) then
c               if (vsf .ge. 0.0) then
c                  fsf = 0.5*(sc1(ijk) + sc1(ijks))
c               else
c                  fsf = sc1(ijk)
c               endif
c            endif
c
c            if (j .eq. nj) then
c               if (vnf .ge. 0.0) then
c                  fnf = sc1(ijk)
c               else
c                  fnf = 0.5*(sc1(ijk) + sc1(ijkn))
c               endif
c            endif
c
c            if (k .eq. 2) then
c               if (wlf .ge. 0.0) then
c                  flf = 0.5*(sc1(ijk) + sc1(ijkl))
c               else
c                  flf = sc1(ijk)
c               endif
c            endif
c
c            if (k .eq. nk) then
c               if (wrf .ge. 0.0) then
c                  frf = sc1(ijk)
c               else
c                  frf = 0.5*(sc1(ijk) + sc1(ijkr))
c               endif
c            endif

            areawf = rdju(i  ,j  )*dkw(k)
            areaef = rdju(i+1,j  )*dkw(k)
            areasf = rdiv(i  ,j  )*dkw(k)
            areanf = rdiv(i  ,j+1)*dkw(k)
            areak  = didjp(i ,j  )
            
            vol = rdidjp(i,j)*dkw(k)

C-----------------------------------------------------------------------
C CONVECTION
C-----------------------------------------------------------------------

C CONVECTION - FLUX TERMS

            gwf = rhowf*uwf*fwf*areawf
            gef = rhoef*uef*fef*areaef
            gsf = rhosf*vsf*fsf*areasf
            gnf = rhonf*vnf*fnf*areanf
            glf = rholf*wlf*flf*areak
            grf = rhorf*wrf*frf*areak

            convf = -(gef-gwf) - (gnf-gsf) - (grf-glf)

C-----------------------------------------------------------------------
C DIFFUSION
C-----------------------------------------------------------------------
                 
C DIFFUSIVE - FLUX TERMS

            galwf = gamwf*areawf
            galef = gamef*areaef
            galsf = gamsf*areasf
            galnf = gamnf*areanf
            gallf = gamlf*areak 
            galrf = gamrf*areak 

            twf = galwf*dfdxwf
            tef = galef*dfdxef
            tsf = galsf*dfdysf
            tnf = galnf*dfdynf
            tlf = gallf*dfdzlf
            trf = galrf*dfdzrf

            difff = (tef-twf) + (tnf-tsf) + (trf-tlf)

C-----------------------------------------------------------------------
C UPDATE ARRAYS
C-----------------------------------------------------------------------

C SCALAR

            h = convf +  difff

            sc1tmp = sc1(ijk)
     &                 + dtime/(vol*rhop)*( 1.5*h - 0.5*hsc1(ijk))

C STORE H(N) FOR NEXT TIME STEP

            hsc1(ijk) = h
!        
!Some amendments to prevent things from happening

               if (sc1tmp .lt. 0.0) sc1tmp = 0.0
               if (sc1tmp .gt. 1.0) sc1tmp = 1.0

!Arbitrary limiters for upper and lower streams
!               if ( iblk .le. 8) sc1tmp = 0.0
               if ( iblk .ge. 1  .and. iblk .le. 8)  sc1tmp = 0.0
               if ( iblk .ge. 25 .and. iblk .le. 32) sc1tmp = 0.0
               if ( iblk .ge. 49 .and. iblk .le. 56) sc1tmp = 1.0
               if ( iblk .ge. 73 .and. iblk .le. 80) sc1tmp = 1.0
               sc1(ijk) = sc1tmp

          end do
        end do
      end do


      return
      end

