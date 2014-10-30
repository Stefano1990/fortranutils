C=======================================================================
C SOURCE TERM FOR PRESSURE EQUATION
C=======================================================================

      subroutine sour_p_block(nip1,njp1,nkp1,io,jo,ko,dtime,
     &                  diu, djv, rdiv, rdju, didjp,
     &                  dkp,
     &                  u,v,w,rho,rhonp1,
     &                  Drhonp1Dt, 
     &                  s,srhot, srhoux,
     &                  gia, gja,
     &                  rdidjp, dkw,blk) 

      implicit none

C ARGUMENTS

      integer      nip1, njp1, nkp1
      integer      io, jo, ko
      real         dtime

      real         diu(nip1,njp1)
      real         djv(nip1,njp1)
      real         rdiv(nip1,njp1)
      real         rdju(nip1,njp1)
      real         didjp(nip1,njp1)

      real         rdidjp(nip1,njp1) 
      real         dkw(nkp1) 

      real         dkp(nkp1)

      real         u(nip1*njp1*nkp1)
      real         v(nip1*njp1*nkp1)
      real         w(nip1*njp1*nkp1)
      real         rho(nip1*njp1*nkp1)
      real         rhonp1(nip1*njp1*nkp1)  
      real         Drhonp1Dt(nip1*njp1*nkp1) 

      real         s(nip1*njp1*nkp1)
      real         srhot(nip1*njp1*nkp1), srhoux(nip1*njp1*nkp1)
      integer      blk  

      real         gia(nkp1), gja(nkp1,nip1)

C LOCAL

      integer      i, j, k
      integer      ni, nj, nk
      integer      joff, ijoff
      integer      ijk
      integer      ijke, ijkn, ijkr

      real         cef, cnf, crf
      real         ddtime
      real         gwf, gef, gsf, gnf, glf, grf

      real         vol   
      real         cdtime 
      real         rdtime  

C COMMON

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------

      ni = nip1 - 1
      nj = njp1 - 1
      nk = nkp1 - 1

C 0.5 IS THE FACTOR IN AVERAGING THE DENSIDY AND IS MOVED TO HERE TO
C REDUCE NUMBER OF OPERATIONS

      ddtime = 2.0/(3.0*dtime) *0.5
      cdtime = 2.0/(3.0*dtime) 
      rdtime = 1.0/dtime  

C-----------------------------------------------------------------------
C DIAGNOSTIC HEADING
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C SET COEFFICIENTS
C-----------------------------------------------------------------------

C INITIALIZE J=2 SURFACE & LOOP THROUGH J

      joff = jo
      do i=2,ni
        cnf = rdiv(i,2)
        ijoff = i + joff
        do k=2,nk
          ijk  = ijoff + ko*(k-1)
cwk          gja(k,i) = cnf*dkp(k)*(rho(ijk)+rho(ijk-jo))*v(ijk)
          gja(k,i) = cnf*dkp(k)*(rhonp1(ijk)+rhonp1(ijk-jo))*v(ijk)
        end do
      end do

      do j=2,nj

C INITIALIZE I = 2 LINE & LOOP THROUGH I

        joff = jo*(j-1)
        cef = rdju(2,j)
        do k=2,nk
          ijk = 2 + joff + ko*(k-1)
cwk          gia(k) = cef*dkp(k)*(rho(ijk)+rho(ijk-1))*u(ijk)
          gia(k) = cef*dkp(k)*(rhonp1(ijk)+rhonp1(ijk-1))*u(ijk)
        end do

        do i=2,ni

C INITIALIZE K = 2 POINT

          cef = rdju(i+1,j)
          cnf = rdiv(i,j+1)
          crf = didjp(i,j)

          ijk  = i + joff + ko
cwk          grf  = crf*(rho(ijk)+rho(ijk-ko))*w(ijk)
          grf  = crf*(rhonp1(ijk)+rhonp1(ijk-ko))*w(ijk)

C LOOP THROUGH K   ?????? SHOULD HAVE A STRIDE OF 1! ???????

          ijoff = i + joff

          do k=2,nk

            ijk  = ijoff + ko*(k-1)
            ijke = ijk + 1
            ijkn = ijk + jo
            ijkr = ijk + ko

            gwf      = gia(k)
cwk            gef      = cef*dkp(k)*(rho(ijk)+rho(ijke))*u(ijke)
            gef      = cef*dkp(k)*(rhonp1(ijk)+rhonp1(ijke))*u(ijke)
            gia(k)   = gef

            gsf      = gja(k,i)
cwk            gnf      = cnf*dkp(k)*(rho(ijk)+rho(ijkn))*v(ijkn)
            gnf      = cnf*dkp(k)*(rhonp1(ijk)+rhonp1(ijkn))*v(ijkn)
            gja(k,i) = gnf

            glf      = grf
cwk            grf      = crf       *(rho(ijk)+rho(ijkr))*w(ijkr)
            grf      = crf       *(rhonp1(ijk)+rhonp1(ijkr))*w(ijkr)

            s(ijk) = -ddtime*(gef - gwf + gnf - gsf + grf - glf)

            srhoux(ijk) =  -ddtime*(gef - gwf + gnf - gsf + grf - glf)
Cwk           Min_flux - Mout_flux
    

            vol = rdidjp(i,j)*dkw(k)
            Drhonp1Dt(ijk) = (rhonp1(ijk)-rho(ijk)) * rdtime        


           s(ijk) = s(ijk) - cdtime*Drhonp1Dt(ijk)*vol  

           srhot(ijk) =  - cdtime*Drhonp1Dt(ijk)*vol
Cwk           Min_flux - Mout_flux - Mchg_rate

          end do

        end do
      end do

      return
      end
