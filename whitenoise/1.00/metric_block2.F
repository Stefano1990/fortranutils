C=======================================================================
C CALCULATE METRIC QUANTITIES
C   dip, diu, div, dix, rdip, rdiv, rdix,
C   djp, dju, djv, djx, rdjp, rdju, rdjx,
C   didjp, didju, rdidju, djdiv, rdjdiv
C FOR A BLOCK
C=======================================================================

      subroutine metric_block(nip1,njp1,nkp1,nip2,njp2,nkp2,
     &  x,y,t,
     &  dip,diu,div,dix,rdip,rdiv,rdix,
     &  djp,dju,djv,djx,rdjp,rdju,rdjx,
     &  hijp,hjip,hkip,hkjp,
     &  hiju,hjiu,hkiu,hkju,
     &  hijv,hjiv,hkiv,hkjv,
     &  rp,ru,rv,
     &  didjp,rdidjp,didju,rdidju,djdiv,rdjdiv,
     &  dkp,dkw,gotrad,iblk)

      implicit none

      integer  nip2, njp2, nkp2
      integer  iblk
      real     x(nip2,njp2), y(nip2,njp2), t(nkp2)
      
      integer  nip1, njp1, nkp1
      real     dip(nip1,njp1), diu(nip1,njp1), div(nip1,njp1),
     &         dix(nip1,njp1), rdip(nip1,njp1), rdiv(nip1,njp1),
     &         rdix(nip1,njp1)
      real     djp(nip1,njp1), dju(nip1,njp1), djv(nip1,njp1),
     &         djx(nip1,njp1), rdjp(nip1,njp1), rdju(nip1,njp1),
     &         rdjx(nip1,njp1)
      real     hijp(nip1,njp1), hjip(nip1,njp1), hkip(nip1,njp1),
     &         hkjp(nip1,njp1)
      real     hiju(nip1,njp1), hjiu(nip1,njp1), hkiu(nip1,njp1),
     &         hkju(nip1,njp1)
      real     hijv(nip1,njp1), hjiv(nip1,njp1), hkiv(nip1,njp1),
     &         hkjv(nip1,njp1)
      real     rp(nip1,njp1), ru(nip1,njp1), rv(nip1,njp1)
      real     didjp(nip1,njp1), rdidjp(nip1,njp1)
      real     didju(nip1,njp1), rdidju(nip1,njp1)
      real     djdiv(nip1,njp1), rdjdiv(nip1,njp1)
      real     dkp(nkp1), dkw(nkp1)
      logical  gotrad

C LOCAL

      integer  i, j, k
      real     x_ij, x_ip1j, x_ijp1, x_ip1jp1, x_ijm1, x_ip1jm1,
     &         x_im1j, x_im1jp1
      real     y_ij, y_ip1j, y_ijp1, y_ip1jp1, y_ijm1, y_ip1jm1,
     &         y_im1j, y_im1jp1
      real     xu_ij, yu_ij, xu_ip1j, yu_ip1j, xu_ijm1, yu_ijm1,
     &         xv_ij, yv_ij, xv_ijp1, yv_ijp1, xv_im1j, yv_im1j,
     &         xp_ij, yp_ij, xp_im1j, yp_im1j, xp_ijm1, yp_ijm1,
     &         ru_ij, rv_ij, rp_ij, rx_ij
      real     rad_on, rad_off, xmax, xmin, ymax, ymin
      real     tiny

      integer  ip1,im1,jp1,jm1,kp1,km1

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      tiny = 1.0e-30

      do j = 1,njp1
        do i = 1,nip1
          dip(i,j)    = 0.
          diu(i,j)    = 0.
          div(i,j)    = 0.
          dix(i,j)    = 0.
          rdip(i,j)   = 0.
          rdiv(i,j)   = 0.
          rdix(i,j)   = 0.
          djp(i,j)    = 0.
          dju(i,j)    = 0.
          djv(i,j)    = 0.
          djx(i,j)    = 0.
          rdjp(i,j)   = 0.
          rdju(i,j)   = 0.
          rdjx(i,j)   = 0.
          hijp(i,j)   = 0.
          hjip(i,j)   = 0.
          hkip(i,j)   = 0.
          hkjp(i,j)   = 0.
          hiju(i,j)   = 0.
          hjiu(i,j)   = 0.
          hkiu(i,j)   = 0.
          hkju(i,j)   = 0.
          hijv(i,j)   = 0.
          hjiv(i,j)   = 0.
          hkiv(i,j)   = 0.
          hkjv(i,j)   = 0.
          rp(i,j)     = 0.
          ru(i,j)     = 0.
          rv(i,j)     = 0.
          didjp(i,j)  = 0.
          rdidjp(i,j) = 0.
          didju(i,j)  = 0.
          rdidju(i,j) = 0.
          djdiv(i,j)  = 0.
          rdjdiv(i,j) = 0.
        end do
      end do

      do k = 1,nkp1
        dkp(k)  = 0.
        dkw(k)  = 0.
      end do

      if (gotrad) then
        rad_on = 1.
      else
        rad_on = 0.
      end if
      rad_off = 1. - rad_on

       xmin = 1.e10
       xmax =-1.e10
       ymin = 1.e10
       ymax =-1.e10
   

C-----------------------------------------------------------------------
C WITHOUT INVOLVING POINTS OUTSIDE OF DOMAIN
C RJ: I'm confused by this - it seemed to do precisely what the
C     comment promises not to do (involve points outside the
C     domain). See original code in metric_block.F, which will
C     generate run time errors if array bounds checking is
C     enabled. 
C     I'm not sure this is really correct now - it depends what
C     the appropriate boundary conditions are for these things.
C-----------------------------------------------------------------------

      do j = 1,njp1
        do i = 1,nip1

          ip1 = min(i+1,nip1)
          im1 = max(i-1,1)
          jp1 = min(j+1,njp1)
          jm1 = max(j-1,1)

          x_ij     = x(i  ,j  )
          x_ip1j   = x(ip1,j  )
          x_ijp1   = x(i  ,jp1)
          x_ip1jp1 = x(ip1,jp1)
          x_im1j   = x(im1,j  )
          x_im1jp1 = x(im1,jp1)
          x_ijm1   = x(i  ,jm1)
          x_ip1jm1 = x(ip1,jm1)
          y_ij     = y(i  ,j  )
          y_ip1j   = y(ip1,j  )
          y_ijp1   = y(i  ,jp1)
          y_ip1jp1 = y(ip1,jp1)
          y_im1j   = y(im1,j  )
          y_im1jp1 = y(im1,jp1)
          y_ijm1   = y(i  ,jm1)
          y_ip1jm1 = y(ip1,jm1)

          xu_ij   = 0.5*(x_ij   + x_ijp1  )
          yu_ij   = 0.5*(y_ij   + y_ijp1  )
          xu_ijm1 = 0.5*(x_ij   + x_ijm1  )
          yu_ijm1 = 0.5*(y_ij   + y_ijm1 )
          xu_ip1j = 0.5*(x_ip1j + x_ip1jp1)
          yu_ip1j = 0.5*(y_ip1j + y_ip1jp1)

          xv_ij   = 0.5*(x_ij   + x_ip1j  )
          yv_ij   = 0.5*(y_ij   + y_ip1j  )
          xv_im1j = 0.5*(x_ij   + x_im1j  )
          yv_im1j = 0.5*(y_ij   + y_im1j  )
          xv_ijp1 = 0.5*(x_ijp1 + x_ip1jp1)
          yv_ijp1 = 0.5*(y_ijp1 + y_ip1jp1)

          xp_ij   = 0.25*(x_ij + x_ip1j + x_ijp1 + x_ip1jp1)
          yp_ij   = 0.25*(y_ij + y_ip1j + y_ijp1 + y_ip1jp1)
          xp_ijm1 = 0.25*(x_ij + x_ip1j + x_ijm1 + x_ip1jm1)
          yp_ijm1 = 0.25*(y_ij + y_ip1j + y_ijm1 + y_ip1jm1)
          xp_im1j = 0.25*(x_ij + x_im1j + x_ijp1 + x_im1jp1)
          yp_im1j = 0.25*(y_ij + y_im1j + y_ijp1 + y_im1jp1)

C CALCULATE ru, rv, rp

          rp_ij = yp_ij*rad_on + rad_off
          ru_ij = yu_ij*rad_on + rad_off
          rv_ij = yv_ij*rad_on + rad_off
          rx_ij = y_ij *rad_on + rad_off

          ru(i,j) = ru_ij
          rv(i,j) = rv_ij
          rp(i,j) = rp_ij

C CALCULATE dip, diu, div, dix, rdip, rdiv, rdix, djp, dju, djv, djx,
C rdjp, rdju, rdjx

          dip(i,j) = sqrt((xu_ip1j-xu_ij)*(xu_ip1j-xu_ij) +
     &                    (yu_ip1j-yu_ij)*(yu_ip1j-yu_ij))
          diu(i,j) = sqrt((xp_ij-xp_im1j)*(xp_ij-xp_im1j) +
     &                    (yp_ij-yp_im1j)*(yp_ij-yp_im1j))
          div(i,j) = sqrt((x_ip1j-x_ij)*(x_ip1j-x_ij) +
     &                    (y_ip1j-y_ij)*(y_ip1j-y_ij))
          dix(i,j) = sqrt((xv_ij-xv_im1j)*(xv_ij-xv_im1j) +
     &                    (yv_ij-yv_im1j)*(yv_ij-yv_im1j))
          rdip(i,j) = rp_ij*dip(i,j)
          rdiv(i,j) = rv_ij*div(i,j)
          rdix(i,j) = rx_ij*dix(i,j)

          djp(i,j) = sqrt((xv_ijp1-xv_ij)*(xv_ijp1-xv_ij) +
     &                    (yv_ijp1-yv_ij)*(yv_ijp1-yv_ij))
          dju(i,j) = sqrt((x_ijp1-x_ij)*(x_ijp1-x_ij) +
     &                    (y_ijp1-y_ij)*(y_ijp1-y_ij))
          djv(i,j) = sqrt((xp_ij-xp_ijm1)*(xp_ij-xp_ijm1) +
     &                    (yp_ij-yp_ijm1)*(yp_ij-yp_ijm1))
          djx(i,j) = sqrt((xu_ij-xu_ijm1)*(xu_ij-xu_ijm1) +
     &                    (yu_ij-yu_ijm1)*(yu_ij-yu_ijm1))
          rdjp(i,j) = rp_ij*djp(i,j)
          rdju(i,j) = ru_ij*dju(i,j)
          rdjx(i,j) = rx_ij*djx(i,j)

c ---- x, y max and min -------

          xmin = min(diu(i,j),xmin)
          xmax = max(diu(i,j),xmax)
          ymin = min(djv(i,j),ymin)
          ymax = max(djv(i,j),ymax)

C CALCULATE didjp, didju, rdidju, djdiv, rdjdiv

          didjp(i,j)  = dip(i,j)*djp(i,j)
          didju(i,j)  = diu(i,j)*dju(i,j)
          djdiv(i,j)  = djv(i,j)*div(i,j)
          rdidjp(i,j) = rp_ij*didjp(i,j)
          rdidju(i,j) = diu(i,j)*rdju(i,j)
          rdjdiv(i,j) = djv(i,j)*rdiv(i,j)

        end do
      end do

c ---

c        write(*,*) 'xmin=', xmin
c        write(*,*) 'xmax=', xmax
c        write(*,*) 'ymin=', ymin
c        write(*,*) 'ymax=', ymax

C CALCULATE dkp, dkw

      do k = 1,nkp1
        kp1 = min(k+1,nkp1)
        km1 = max(k-1,1)
        dkp(k)  = t(kp1) - t(k)
        dkw(k) = (t(kp1) - t(km1))*0.5
      end do

C CALCULATE hijp, hjip, hkip, hkjp; hiju, hjiu, hkiu, hkju;
C           hijv, hjiv, hkiv, hkjv

      do j = 2,njp1-1
        do i = 2,nip1-1

          hijp(i,j) = (div(i,j+1)-div(i,j))/(didjp(i,j)+tiny)
          hjip(i,j) = (dju(i+1,j)-dju(i,j))/(didjp(i,j)+tiny)
          hkip(i,j) = (ru(i+1,j) -ru(i,j) )/(rdip(i,j)+tiny)
          hkjp(i,j) = (rv(i,j+1) -rv(i,j) )/(rdjp(i,j)+tiny)

          hiju(i,j) = (dix(i,j+1)-dix(i,j))/(didju(i,j)+tiny)
          hjiu(i,j) = (djp(i,j)-djp(i-1,j))/(didju(i,j)+tiny)
          hkiu(i,j) = (rp(i,j) -rp(i-1,j) )/(ru(i,j)*diu(i,j)+tiny)
          hkju(i,j) = (y_ijp1 - y_ij)*rad_on/(rdju(i,j)+tiny)

          hijv(i,j) = (dip(i,j)-dip(i,j-1))/(djdiv(i,j)+tiny)
          hjiv(i,j) = (djx(i+1,j)-djx(i,j))/(djdiv(i,j)+tiny)
          hkiv(i,j) = (y_ip1j - y_ij)*rad_on/(rdiv(i,j)+tiny)
          hkjv(i,j) = (rp(i,j) - rp(i,j-1))/(rv(i,j)*djv(i,j)+tiny)

        end do
      end do

      return
      end
