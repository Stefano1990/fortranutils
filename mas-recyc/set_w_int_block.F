C=======================================================================
C SET INTERMEDIATE W
C=======================================================================

C EQUATIONS FORMULATED WITH RADIUS REMOVED FROM DENOMINATOR

C SET W*(N+1), HW(N) FROM HW(N-1), U(N), V(N), W(N) & P(N-1)

      subroutine set_w_int_block(nip1,njp1,nkp1,io,jo,ko,dtime,
     &                     dip, diu, rdiv,
     &                     djp, djv, rdju,
     &                     didjp, rdidjp,
     &                     hkip, hkjp, hkiu, hkjv,
     &                     rp, ru,
     &                     dkp, dkw,
     &                     w_int,hw,u,v,w,p,rho,rhonp1,gam, 
     &                     tw_sfk, tw_nfk, tw_wfk, tw_efk,
     &                     wallwf,wallef,wallsf,wallnf,walllf,wallrf,
     &                     kstart,iststep,visc,ylaminar,kappa,cwfun)

C PARAMETERS

C ARGUMENTS

      integer      nip1, njp1, nkp1
      integer      io, jo, ko
      integer      kstart
      real         dtime

      real         dip(nip1,njp1), diu(nip1,njp1), rdiv(nip1,njp1)
      real         djp(nip1,njp1), djv(nip1,njp1), rdju(nip1,njp1)
      real         didjp(nip1,njp1), rdidjp(nip1,njp1)
      real         hkip(nip1,njp1), hkjp(nip1,njp1)
      real         hkiu(nip1,njp1), hkjv(nip1,njp1)
      real         rp(nip1,njp1), ru(nip1,njp1)
      real         dkp(nkp1), dkw(nkp1)

      real         w_int(nip1*njp1*nkp1)
      real         hw(nip1*njp1*nkp1)
      real         u(nip1*njp1*nkp1)
      real         v(nip1*njp1*nkp1)
      real         w(nip1*njp1*nkp1)
      real         p(nip1*njp1*nkp1)
      real         rho(nip1*njp1*nkp1)
      real         rhonp1(nip1*njp1*nkp1)  
      real         gam(nip1*njp1*nkp1)

      real         tw_sfk(nkp1,nip1), tw_nfk(nkp1,nip1)
      real         tw_wfk(nkp1,njp1), tw_efk(nkp1,njp1)

      logical      wallwf, wallef, wallsf, wallnf, walllf, wallrf
      logical      turbul

      integer      iststep
      real         visc
      real         ylaminar, kappa, cwfun

C LOCAL

      integer      i, j, k
      integer      ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr
      integer      ijksw, ijknw, ijkwl, ijkwr
      integer      ijkse, ijksl, ijksr
      integer      ijkel, ijknl
      integer      ni, nj, nk
      
      real         rhowf, rhoef, rhosf, rhonf, rholf, rhorf, rhow
      real         rhonp1w  
      real         gamwf, gamef, gamsf, gamnf, gamlf, gamrf
      real         galwf, galef, galsf, galnf, gallf, galrf
      real         areawf, areaef, areasf, areanf, areak
      real         gwf, gef, gsf, gnf, glf, grf
      real         cwf, cef, csf, cnf, clf, crf
      real         twf, tef, tsf, tnf, tlf, trf
      real         uwf, uef, ulf, urf
      real         vsf, vnf, vlf, vrf
      real         wlf, wrf, wwf, wef,  wsf, wnf
      real         tikwf, tikef, tjksf, tjknf, tkklf, tkkrf
      real         h
      real         st1, st2, stden
      real         vol
      real         conf, conc, diffc, difff, diffc1, diffc2
      real         dist, tw, tw_mean, ut_mean, yplus_mean
      real         uloglaw
      real         di, dj
      real         radp, radw, rade, rads, radn
      real         radwf, radef, radsf, radnf
      real         dwdy, dvdz, dwdx, dudz

C COMMON

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------

      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C GATHERING STATISTICS

      st1   = iststep - 1
      st2   = iststep - 2
      stden = 1.0/iststep
      
C-----------------------------------------------------------------------
C SET W* & HW
C-----------------------------------------------------------------------

      do j=2,nj
        do i=2,ni

C SET I,J METRICS

          di = dip(i,j)
          dj = djp(i,j)

          radp = rp(i  ,j  )
          radw = rp(i-1,j  )
          rade = rp(i+1,j  )
          rads = rp(i  ,j-1)
          radn = rp(i  ,j+1)

          do k=kstart,nk

            ijk  = i + jo*(j-1) + ko*(k-1)

            ijkw = ijk - io
            ijke = ijk + io
            ijks = ijk - jo
            ijkn = ijk + jo
            ijkl = ijk - ko
            ijkr = ijk + ko

            ijksw = ijks - io
            ijknw = ijkn - io
            ijkwl = ijkw - ko
            ijkwr = ijkw + ko

            ijkse = ijks + io
            ijksl = ijks - ko
            ijksr = ijks + ko

            ijkel = ijke - ko
            ijknl = ijkn - ko
            ijkne = ijkn + io

C-----------------------------------------------------------------------
C FACE AND CELL AVERAGE VALUES
C-----------------------------------------------------------------------

            rhowf = 0.25*(rho(ijk)+rho(ijkl)+rho(ijkw)+rho(ijkwl))
            rhoef = 0.25*(rho(ijk)+rho(ijkl)+rho(ijke)+rho(ijkel))
            rhosf = 0.25*(rho(ijk)+rho(ijkl)+rho(ijks)+rho(ijksl))
            rhonf = 0.25*(rho(ijk)+rho(ijkl)+rho(ijkn)+rho(ijknl))
            rholf = rho(ijkl)
            rhorf = rho(ijk )
            rhow  = 0.5*(rho(ijk)+rho(ijkl))

            rhonp1w = 0.5*(rhonp1(ijk)+rhonp1(ijkl))


            gamwf = visc
     &            + .25*(gam(ijk)+gam(ijkl)+gam(ijkw)+gam(ijkwl))*rhowf
            gamef = visc
     &            + .25*(gam(ijk)+gam(ijkl)+gam(ijke)+gam(ijkel))*rhoef
            gamsf = visc
     &            + .25*(gam(ijk)+gam(ijkl)+gam(ijks)+gam(ijksl))*rhosf
            gamnf = visc
     &            + .25*(gam(ijk)+gam(ijkl)+gam(ijkn)+gam(ijknl))*rhonf
            gamlf = visc + gam(ijkl)*rholf
            gamrf = visc + gam(ijk )*rhorf

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

C CURVATURE TERMS INCORPORATED INTO COEFFICIENTS - NON-CONSERVATIVE!

            uwf = 0.5*(u(ijk ) + u(ijkl ))
            uef = 0.5*(u(ijke) + u(ijkel))

            vsf = 0.5*(v(ijk ) + v(ijkl ))
            vnf = 0.5*(v(ijkn) + v(ijknl))

            wlf = 0.5*(w(ijk ) + w(ijkl ))
            wrf = 0.5*(w(ijkr) + w(ijk  ))

            gwf = 0.5*(w(ijk)+w(ijkw))*uwf*rhowf*areawf
            gef = 0.5*(w(ijk)+w(ijke))*uef*rhoef*areaef
            
            gsf = 0.5*(w(ijk)+w(ijks))*vsf*rhosf*areasf
            gnf = 0.5*(w(ijk)+w(ijkn))*vnf*rhonf*areanf
            
            glf = wlf*wlf*rholf*areak
            grf = wrf*wrf*rhorf*areak
                        
            conf = -(gef-gwf) - (gnf-gsf) - (grf-glf)

C CURVATURE TERMS 

	    u_wc = 0.25*(u(ijk)+u(ijke)+u(ijkl)+u(ijkel))
	    v_wc = 0.25*(v(ijk)+v(ijkn)+v(ijkl)+v(ijknl))
	    
	    conc = -vol*rhow*(w(ijk)*(u_wc*hkip(i,j)+v_wc*hkjp(i,j)))

c	write(*,*) 'I am in the 3 of set_w_int_block'


C-----------------------------------------------------------------------
C DIFFUSION
C-----------------------------------------------------------------------
                 
C DIFFUSION - FLUX TERMS

C CURVATURE TERMS INCORPORATED INTO COEFFICIENTS - NON-CONSERVATIVE!

C - PACE INTERPOLATES (GAM*RAD**2) FOR A FACE. WHY???????

C????????????????????????
C????????????????????????
C - ABSORBING C.C. TERMS => "FLUX TERMS" ARE NOT STRESSES! 
C????????????????????????
C????????????????????????

            galwf = gamwf*areawf/diu(i  ,j  )
            galef = gamef*areaef/diu(i+1,j  )
            galsf = gamsf*areasf/djv(i  ,j  )
            galnf = gamnf*areanf/djv(i  ,j+1)
            gallf = gamlf*areak /(radp*dkp(k-1))
            galrf = gamrf*areak /(radp*dkp(k  ))

c	write(*,*) 'I am in the 4 of set_w_int_block'
            radwf = 0.5*(radp + radw)
            radef = 0.5*(radp + rade)
            radsf = 0.5*(radp + rads)
            radnf = 0.5*(radp + radn)

            ulf = 0.5*(u(ijkl) + u(ijkel))
            urf = 0.5*(u(ijk ) + u(ijke ))

            vlf = 0.5*(v(ijkl) + v(ijknl))
            vrf = 0.5*(v(ijk ) + v(ijkn ))

            wwf = 0.5*(w(ijk ) + w(ijkw ))
            wef = 0.5*(w(ijk ) + w(ijke ))

            wsf = 0.5*(w(ijk ) + w(ijks ))
            wnf = 0.5*(w(ijk ) + w(ijkn ))

            twf = galwf*radwf*(w(ijk )/radp - w(ijkw)/radw)
            tef = galef*radef*(w(ijke)/rade - w(ijk )/radp)

            tsf = galsf*(w(ijk ) - w(ijks))
            tnf = galnf*(w(ijkn) - w(ijk ))

            tlf = 2.0*gallf*(w(ijk )-w(ijkl))
            trf = 2.0*galrf*(w(ijkr)-w(ijk ))

c	write(*,*) 'I am in the 5 of set_w_int_block'
            cwf = gamwf*areawf*((u(ijk )-u(ijkl ))/dkw(k)/radwf
     &            -wwf*hkiu(i,j))
            cef = gamef*areaef*((u(ijke)-u(ijkel))/dkw(k)/radef
     &            -wef*hkiu(i+1,j))

c	write(*,*) 'I am in the 6 of set_w_int_block'

          if(j.eq.2) then
           csf = 0.0
          else
           csf = gamsf*areasf*((v(ijk )-v(ijkl ))/dkw(k)/radsf
     &            -wsf*hkjv(i,j))
          endif
           cnf = gamnf*areanf*((v(ijkn)-v(ijknl))/dkw(k)/radnf
     &            -wnf*hkjv(i,j+1))
c	write(*,*) 'I am in the 7 of set_w_int_block'

            clf = 2.0*gamlf*areak*(ulf*hkip(i,j) + vlf*hkjp(i,j))
            crf = 2.0*gamrf*areak*(urf*hkip(i,j) + vrf*hkjp(i,j))

c I CHANGED THIS
c           clf = 2.0*gamlf*areak*radp*(ulf*hkip(i,j) + vlf*hkjp(i,j))
c           crf = 2.0*gamrf*areak*radp*(urf*hkip(i,j) + vrf*hkjp(i,j))

c --- Curvature terms of diffusion -----

	    dwdy = (wnf-wsf)/djp(i,j)
       	    dvdz = (vrf-vlf)/(radp*dkw(k))
            dwdx = (wef-wwf)/dip(i,j)
            dudz = (urf-ulf)/(radp*dkw(k))

            diffc1= hkip(i,j)*(dwdx+dudz - w(ijk)*hkip(i,j))*vol*gamrf
            diffc2= hkjp(i,j)*(dwdy+dvdz - w(ijk)*hkjp(i,j))*vol*gamrf
            diffc = diffc1 + diffc2

C DIFFUSIVE STRESSES ON FACES - INCORRECT (CC TERMS!) ???????

            tikwf = twf + cwf
            tikef = tef + cef
            tjksf = tsf + csf
            tjknf = tnf + cnf
            tkklf = tlf + clf
            tkkrf = trf + crf

C-----------------------------------------------------------------------
C MODIFICATION FOR WALL: OVERWRITE FACE FLUXES & UPDATE MEAN VALUE
C-----------------------------------------------------------------------

C  - USE STATISTICAL MEAN FOR U????  
C  - PROBLEM WITH SWITCH AT YLAMINAR!!!    
C  - STORE UT AND NOT TW
C  - TWALL HAS RHO MISSING
C  - NO TREATMENT FOR WF & EF (EXCEPT SKIPPING I=2)
C  - COMMENT NOT CODED FOR MOVING WALL ----- ADVANCED BY GEFENG TANG

c -- tangtangtang ---------
c	turbul = .true.
c	if(turbul) then
c -- tangtangtang --------

C WEST WALL

            if (wallwf .and. i.eq.2) then

              dist       = 0.5*dip(i,j)
              tw_mean    = abs(tw_wfk(k,j))
              ut_mean    = sqrt(tw_mean/rhowf)
              yplus_mean = dist*ut_mean*rhowf/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = w(ijk)*tw_mean/uloglaw
              else
                tw      = visc*w(ijk)/dist
              end if

              tikwf       =  areawf*tw
c              tw_wfk(k,j) =  (st1*(tw+tw_mean) - st2*tw)*stden
              tw_wfk(k,j) =  (st1*(tw+tw_wfk(k,j)) - st2*tw)*stden

            end if

C EAST WALL

            if (wallef .and. i.eq.ni) then

              dist       = 0.5*dip(i,j)
              tw_mean    = abs(tw_efk(k,j))
              ut_mean    = sqrt(tw_mean/rhoef)
              yplus_mean = dist*ut_mean*rhoef/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = w(ijk)*tw_mean/uloglaw
              else
                tw      = visc*w(ijk)/dist
              end if

              tikef       = -areaef*tw
c              tw_efk(k,j) = (st1*(-tw+tw_mean) - st2*(-tw))*stden
              tw_efk(k,j) = (st1*(-tw+tw_efk(k,j)) - st2*(-tw))*stden

            end if

C SOUTH WALL

c-- for rough wall ----

c           cwfun_in_b = 0.50
c           cwfun_in_t = 1.90
c

            if (wallsf .and. j.eq.2) then

              dist       = 0.5*djp(i,j)
              tw_mean    = abs(tw_sfk(k,i))
              ut_mean    = sqrt(tw_mean/rhosf)
              yplus_mean = dist*ut_mean*rhosf/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = w(ijk)*tw_mean/uloglaw
              else
c --- for moving wall ----

                w_wall = 0.0
                tw      = visc*(w(ijk)-w_wall)/dist
              end if

              tjksf       =  areasf*(tw-w_wall/0.005*visc)
c              tw_sfk(k,i) =  (st1*(tw+tw_mean) - st2*tw)*stden
              tw_sfk(k,i) =  (st1*(tw+tw_sfk(k,i)) - st2*tw)*stden

            end if

C NORTH WALL

c	   write(*,*) 'moving wall ---- I am here'

            if (wallnf .and. j.eq.nj) then

              dist       = 0.5*djp(i,j)
              tw_mean    = abs(tw_nfk(k,i))
              ut_mean    = sqrt(tw_mean/rhonf)
              yplus_mean = dist*ut_mean*rhonf/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = w(ijk)*tw_mean/uloglaw
              else
c --- for moving wall ----

                w_wall = 0.0
                tw      = visc*(w(ijk)-w_wall)/dist
              end if

              tjknf       = -areanf*(tw+w_wall/0.015*visc)
c              tw_nfk(k,i) = (st1*(-tw+tw_mean) - st2*(-tw))*stden
              tw_nfk(k,i) = (st1*(-tw+tw_nfk(k,i)) - st2*(-tw))*stden

            end if

c	endif

C-----------------------------------------------------------------------
C UPDATE ARRAYS
C-----------------------------------------------------------------------

C INTERMEDIATE VELOCITY

            difff = tikef - tikwf + tjknf - tjksf + tkkrf- tkklf
            h = conf + conc +  difff + diffc


            w_int(ijk) = w(ijk)*(rhow/rhonp1w)
     &                 + dtime/(vol*rhonp1w)*( 1.5*h - 0.5*hw(ijk) +
     &                             0.5*areak*(p(ijk) - p(ijkl)) )


C STORE H(N) FOR NEXT TIME STEP

            hw(ijk) = h

          end do
        end do
      end do



      return
      end

