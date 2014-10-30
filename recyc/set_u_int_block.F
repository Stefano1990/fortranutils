C=======================================================================
C SET INTERMEDIATE U
C=======================================================================

C SET U*(N+1), HU(N) FROM HU(N-1), U(N), V(N), W(N) & P(N-1)

      subroutine set_u_int_block(nip1,njp1,nkp1,io,jo,ko,io_ij,jo_ij,
     &                     dtime, dip, diu, dix, rdix,
     &                     djp, dju, djv, djx, rdjp, rdju,
     &                     didjp, didju, rdidju,
     &                     hijp, hiju, hjiu, hkiu, hkju,
     &                     rp, ru,
     &                     dkp, dkw,
     &                     u_int,hu,u,v,w,p,rho,rhonp1,gam,
     &                     tw_sfi, tw_nfi, tw_lfi, tw_rfi,
     &                     wallwf,wallef,wallsf,wallnf,walllf,wallrf,
     &                     istart,iststep,visc,ylaminar,kappa,cwfun,iva)

C PARAMETERS

C ARGUMENTS

      integer      nip1, njp1, nkp1
      integer      io, jo, ko, io_ij, jo_ij
      real         dtime
      integer      istart, iva

      real         dip(nip1*njp1), diu(nip1*njp1), dix(nip1*njp1)
      real         rdix(nip1*njp1)
      real         djp(nip1*njp1), dju(nip1*njp1), djv(nip1*njp1)
      real         djx(nip1*njp1)
      real         rdjp(nip1*njp1), rdju(nip1*njp1)
      real         didjp(nip1*njp1), didju(nip1*njp1), rdidju(nip1*njp1)
      real         hijp(nip1*njp1), hiju(nip1*njp1), hjiu(nip1*njp1)
      real         hkiu(nip1*njp1), hkju(nip1*njp1)
      real         rp(nip1*njp1), ru(nip1*njp1)

      real         dkp(nkp1), dkw(nkp1)

      real         u_int(nip1*njp1*nkp1)
      real         hu(nip1*njp1*nkp1)
      real         u(nip1*njp1*nkp1)
      real         v(nip1*njp1*nkp1)
      real         w(nip1*njp1*nkp1)
      real         p(nip1*njp1*nkp1)
      real         rho(nip1*njp1*nkp1)
      real         rhonp1(nip1*njp1*nkp1) 
      real         gam(nip1*njp1*nkp1)

      real         tw_sfi(nkp1*nip1), tw_nfi(nkp1*nip1)
      real         tw_lfi(nip1*njp1), tw_rfi(nip1*njp1)

      logical      wallwf, wallef, wallsf, wallnf, walllf, wallrf

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
      integer      joff_ij, joff, ijoff, ioff_ki
      integer      ij, ijw, ijs, ijn, ijnw
      integer      ki

      real         rhowf, rhoef, rhosf, rhonf, rholf, rhorf, rhou
      real         rhonp1u 
      real         gamwf, gamef, gamsf, gamnf, gamlf, gamrf
      real         galwf, galef, galsf, galnf, gallf, galrf
      real         areawf, areaef, areasf, areanf, areak
      real         gwf, gef, gsf, gnf, glf, grf
      real         cwf, cef, csf, cnf, clf, crf
      real         twf, tef, tsf, tnf, tlf, trf
      real         uwf, uef
      real         vwf, vef, vsf, vnf
      real         wlf, wrf
      real         tiiwf, tiief, tjisf, tjinf, tkilf, tkirf
      real         udisf, udinf
      real         vdjwf, vdjef
      real         h
      real         st1, st2, stden
      real         vol, rdkp, rdku
      real         conf, conc
      real         cc1, cc2, cc3
      real         diffc
      real         dc1, dc2, dc3
      real         dist, tw, tw_mean, ut_mean, yplus_mean
      real         uloglaw
      real         va, wa, rhoa, gama
      real         di, dj

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
C SET U* & HU
C-----------------------------------------------------------------------

      do j=2,nj

        joff    = 1 + jo   *(j-1)
        joff_ij = 1 + jo_ij*(j-1)

        do i=istart,ni

          ijoff = io*(i-1) + joff

          ij   = io_ij*(i-1) + joff_ij
          ijw  = ij  - io_ij
          ijs  = ij  - jo_ij
          ijn  = ij  + jo_ij
          ijnw = ijw + jo_ij

          ioff_ki = nkp1*(i-1)

C SET I,J METRICS

          di = diu(ij)
          dj = dju(ij)

          do k=2,nk

            ijk = ko*(k-1) + ijoff
            ki  = k        + ioff_ki

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

            rdkp = rp(ij)*dkp(k)
            rdku = ru(ij)*dkp(k)

C-----------------------------------------------------------------------
C FACE AND CELL AVERAGE VALUES
C-----------------------------------------------------------------------

            va   = 0.25*(v(ijk) + v(ijkn) + v(ijkw) + v(ijknw))
            wa   = 0.25*(w(ijk) + w(ijkr) + w(ijkw) + w(ijkwr))
            rhoa = 0.5 *(rho(ijkw) + rho(ijk))
            gama = 0.5 *(gam(ijkw) + gam(ijk))

C THE FOLLOWING CALCULATION OF FACE VALUES OF rho REQUIRES VIRTUAL
C VALUES OF rho NEAR BOUNDARIES IF NO IF-THEN STATEMENTS ARE USED.
C HOWEVER, FOR CONSTANT rho, THE FOLLOWING FORMULA ARE CORRENT PROVIDED
C THAT THE POINTS FOR rho OUTSIDE OF BOUNDARIES ARE ASSINGED THE VALUE.
C IF rho IS NOT CONSTANT, set_virt_bc_block NEEDS TO BE MADE MORE
C GENERAL TO SUPPLY THE VIRTUAL VALUES.

            rhowf = rho(ijkw)
            rhoef = rho(ijk )
            rhosf = 0.25*(rho(ijk)+rho(ijkw)+rho(ijks)+rho(ijksw))
            rhonf = 0.25*(rho(ijk)+rho(ijkw)+rho(ijkn)+rho(ijknw))
            rholf = 0.25*(rho(ijk)+rho(ijkw)+rho(ijkl)+rho(ijkwl))
            rhorf = 0.25*(rho(ijk)+rho(ijkw)+rho(ijkr)+rho(ijkwr))
            rhou  = 0.5 *(rho(ijk)+rho(ijkw))
            rhonp1u = 0.5 *(rhonp1(ijk)+rhonp1(ijkw))



C NOTE: gam MAY NOT BE DEFINED. FOR INSTANCE, ON THE SOUTH WALL, gam AT
C ijks AND ijksw ARE NOT DEFINED, AND SO gamsf CALCULATED IS USELESS.
C HOWEVER, THIS DOES NOT MATTER, SINCE THE QUANTITY, tjisf, CALCULATED
C BY USING gamsf WILL BE OVERWRITTEN LATER.

            gamwf = visc + gam(ijkw)*rhowf
            gamef = visc + gam(ijk )*rhoef
            gamsf = visc
     &            + .25*(gam(ijk)+gam(ijkw)+gam(ijks)+gam(ijksw))*rhosf
            gamnf = visc
     &            + .25*(gam(ijk)+gam(ijkw)+gam(ijkn)+gam(ijknw))*rhonf
            gamlf = visc
     &            + .25*(gam(ijk)+gam(ijkw)+gam(ijkl)+gam(ijkwl))*rholf
            gamrf = visc
     &            + .25*(gam(ijk)+gam(ijkw)+gam(ijkr)+gam(ijkwr))*rhorf

            areawf = rdjp(ijw)*dkp(k)
            areaef = rdjp(ij )*dkp(k)
            areasf = rdix(ij )*dkp(k)
            areanf = rdix(ijn)*dkp(k)
            areak  = didju(ij)

            vol = rdidju(ij)*dkp(k)

C-----------------------------------------------------------------------
C CONVECTION
C-----------------------------------------------------------------------

C CONVECTION - FLUX TERMS

            uwf = 0.5*(u(ijk) + u(ijkw))
            uef = 0.5*(u(ijk) + u(ijke))

            vsf = 0.5*(v(ijk ) + v(ijkw ))
            vnf = 0.5*(v(ijkn) + v(ijknw))

            wlf = 0.5*(w(ijk ) + w(ijkw ))
            wrf = 0.5*(w(ijkr) + w(ijkwr))

            gwf =                  uwf*uwf*rhowf*areawf
            gef =                  uef*uef*rhoef*areaef

            gsf = 0.5*(u(ijk)+u(ijks))*vsf*rhosf*areasf
            gnf = 0.5*(u(ijk)+u(ijkn))*vnf*rhonf*areanf
            glf = 0.5*(u(ijk)+u(ijkl))*wlf*rholf*areak
            grf = 0.5*(u(ijk)+u(ijkr))*wrf*rhorf*areak

            conf = -(gef-gwf) - (gnf-gsf) - (grf-glf)

C CONVECTION - CURVATURE TERMS

            cc1 =  hjiu(ij)*rhoa*(va**2)  *vol
            cc2 =  hkiu(ij)*rhoa*(wa**2)  *vol
            cc3 = -hiju(ij)*rhoa*va*u(ijk)*vol

            conc = cc1 + cc2 + cc3  

C-----------------------------------------------------------------------
C DIFFUSION
C-----------------------------------------------------------------------

C DIFFUSION - FLUX TERMS

C?????????????????? SORT OUT!!! ???????????
C?????????????????? SORT OUT!!! ???????????
C NOTE: MY EQUATIONS HAVE 2.0*DE1 BUT PACE HAS 1.0*DE1 !!!
C?????????????????? SORT OUT!!! ???????????
C?????????????????? SORT OUT!!! ???????????

            galwf = gamwf*areawf/dip(ijw)
            galef = gamef*areaef/dip(ij )
            galsf = gamsf*areasf/djx(ij )
            galnf = gamnf*areanf/djx(ijn)
            gallf = gamlf*areak /(rp(ij)*dkw(k  ))
            galrf = gamrf*areak /(rp(ij)*dkw(k+1))

            vwf = 0.5*(v(ijkw) + v(ijknw))
            vef = 0.5*(v(ijk ) + v(ijkn ))

            twf = galwf*2.0*(u(ijk )-u(ijkw))
            tef = galef*2.0*(u(ijke)-u(ijk ))

C NOTE: diu(i,j-1) IS OUTSIDE OF BOUNDARY FOR j=2 AND THEREFORE IS NOT
C DEFINED. THE VALUE OF diu FOR THIS CASE IS TAKEN TO BE THAT ON THE
C WALL. THE SAME APPLIES TO diu(i,j+1) FOR j=nj. HOWEVER, ON WALLS, 
C tjisf AND tjinf WILL BE OVERWRITTEN ANYWAY, SO ACTUALLY WHETHER diu
C IS DEFINED OR NOT DOES NOT MATTER, AS LONG AS IT'S NOT ZERO.

            tsf = galsf*dix(ij )*(u(ijk )/diu(ij ) -
     &                            u(ijks)/diu(ijs)   )
            tnf = galnf*dix(ijn)*(u(ijkn)/diu(ijn) -
     &                            u(ijk )/diu(ij )   )

            tlf = gallf*(u(ijk )-u(ijkl))
            trf = galrf*(u(ijkr)-u(ijk ))

            cwf = 2.0*gamwf*areawf*vwf*hijp(ijw)
            cef = 2.0*gamef*areaef*vef*hijp(ij )

            csf = gamsf*areasf*djx(ij )*(v(ijk  )/djv(ij ) -
     &                                   v(ijkw )/djv(ijw )   ) /
     &                                   dix(ij )
            cnf = gamnf*areanf*djx(ijn)*(v(ijkn )/djv(ijn) -
     &                                   v(ijknw)/djv(ijnw)   ) /
     &                                   dix(ijn)

            clf = gamlf*areak*ru(ij)*(w(ijk  )/rp(ij ) - 
     &                                w(ijkw )/rp(ijw)   ) / diu(ij)
            crf = gamrf*areak*ru(ij)*(w(ijkr )/rp(ij ) - 
     &                                w(ijkwr)/rp(ijw)   ) / diu(ij)

C DIFFUSIVE STRESSES ON FACES

            tiiwf = twf + cwf
            tiief = tef + cef
            tjisf = tsf + csf
            tjinf = tnf + cnf
            tkilf = tlf + clf
            tkirf = trf + crf

C DIFFUSION - CURVATURE TERMS

            vsf = 0.5*(v(ijk ) + v(ijkw ))
            vnf = 0.5*(v(ijkn) + v(ijknw))

            dc1 = -2.0*gama*vol*hjiu(ij)*((vnf-vsf)/dju(ij) + 
     &                                      u(ijk)*hjiu(ij)     )

            wlf = 0.5*(w(ijk ) + w(ijkw ))
            wrf = 0.5*(w(ijkr) + w(ijkwr))

            dc2 = -2.0*gama*vol*hkiu(ij)*(u(ijk)*hkiu(ij) +
     &                                    va    *hkju(ij) +
     &                                     (wrf-wlf)/rdkp)

            udisf = 0.5*(u(ijk)+u(ijks))/dix(ij )
            udinf = 0.5*(u(ijk)+u(ijkn))/dix(ijn)

            vdjwf = 0.5*(v(ijkw)+v(ijknw))/djp(ijw)
            vdjef = 0.5*(v(ijk )+v(ijkn ))/djp(ij )

            dc3 = gama*vol*hiju(ij)*(di*(udinf-udisf)/dj +
     &                               dj*(vdjef-vdjwf)/di   )

            diffc = dc1 + dc2 + dc3

C-----------------------------------------------------------------------
C MODIFICATION FOR WALL: OVERWRITE FACE FLUXES & UPDATE MEAN VALUE
C-----------------------------------------------------------------------

C  - USE STATISTICAL MEAN FOR U????  
C  - PROBLEM WITH SWITCH AT YLAMINAR!!!    
C  - STORE UT AND NOT TW
C  - TWALL HAS RHO MISSING
C  - NO TREATMENT FOR WF & EF (EXCEPT SKIPPING I=2)



C SOUTH WALL

            if (wallsf .and. j.eq.2) then

              dist       = 0.5*dju(ij)
              tw_mean    = abs(tw_sfi(ki))
              ut_mean    = sqrt(tw_mean/rhosf)
c             write(*,*) 'ut_mean j=2 is**********',ut_mean
              yplus_mean = dist*ut_mean*rhosf/visc


              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = u(ijk)*tw_mean/uloglaw
              else
                tw      = visc*u(ijk)/dist
              end if

              tjisf      =  areasf*tw
c              tw_sfi(ki) =  (st1*(tw+tw_mean) - st2*tw)*stden
              tw_sfi(ki) =  (st1*(tw+tw_sfi(ki)) - st2*tw)*stden

            end if

C NORTH WALL

            if (wallnf .and. j.eq.nj) then

              dist       = 0.5*dju(ij)
              tw_mean    = abs(tw_nfi(ki))
              ut_mean    = sqrt(tw_mean/rhonf)
              yplus_mean = dist*ut_mean*rhonf/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = u(ijk)*tw_mean/uloglaw
              else
                tw      = visc*u(ijk)/dist
              end if

              tjinf      = -areanf*tw
c              tw_nfi(ki) = (st1*(-tw+tw_mean) - st2*(-tw))*stden
              tw_nfi(ki) = (st1*(-tw+tw_nfi(ki)) - st2*(-tw))*stden

            end if

C LEFT WALL

            if (walllf .and. k.eq.2) then

              dist       = 0.5*rdku
              tw_mean    = abs(tw_lfi(ij))
              ut_mean    = sqrt(tw_mean/rholf)
              yplus_mean = dist*ut_mean*rholf/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = u(ijk)*tw_mean/uloglaw
              else
                tw      = visc*u(ijk)/dist
              end if

              tkilf      =  areak*tw
c              tw_lfi(ij) =  (st1*(tw+tw_mean) - st2*tw)*stden
              tw_lfi(ij) =  (st1*(tw+tw_lfi(ij)) - st2*tw)*stden

            end if

C RIGHT WALL

            if (wallrf .and. k.eq.nk) then

              dist       = 0.5*rdku
              tw_mean    = abs(tw_rfi(ij))
              ut_mean    = sqrt(tw_mean/rhorf)
              yplus_mean = dist*ut_mean*rhorf/visc

              if (yplus_mean.gt.ylaminar) then
                uloglaw = ut_mean/kappa*(log(yplus_mean) + cwfun)
                tw      = u(ijk)*tw_mean/uloglaw
              else
                tw      = visc*u(ijk)/dist
              end if

              tkirf      = -areak*tw 
c              tw_rfi(ij) = (st1*(-tw+tw_mean) - st2*(-tw))*stden
              tw_rfi(ij) = (st1*(-tw+tw_rfi(ij)) - st2*(-tw))*stden

            end if

C-----------------------------------------------------------------------
C UPDATE ARRAYS
C-----------------------------------------------------------------------

C INTERMEDIATE VELOCITY

            h = conf + conc + diffc +
     &          tiief - tiiwf + tjinf - tjisf + tkirf - tkilf


            u_int(ijk) = u(ijk)*(rhou/rhonp1u)
     &                 + dtime/(rhonp1u*vol)*( 1.5*h - 0.5*hu(ijk) +
     &                     0.5*(p(ijk) - p(ijkw))*vol/diu(ij) )


C STORE H(N) FOR NEXT TIME STEP

            hu(ijk) = h

          end do
        end do
      end do

      return
      end

