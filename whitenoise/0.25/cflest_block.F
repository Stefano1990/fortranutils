C=======================================================================
C CALCULATE THE LARGEST CFL AND DFS NUMBERS IN A BLOCK
C=======================================================================

      subroutine cflest_block(u,v,w,sc1,
     &  rho,Temperature, 
     &  dtime,visc,dip,djp,rp,dkw,gam,
     &  istr,jstr,kstr,nip1,njp1,nkp1,cflmax,dfsmax,ijkmaxcfl,
     &  ijkmaxdfs,uvw_min,uvw_max,sc1_min,sc1_max,
     &  rho_min,rho_max, T_min, T_max, 
     &  ijkmaxsc1,ijkminsc1,
     &  ijkmaxrho,ijkminrho,ijkmaxT,ijkminT) 

      implicit none

#include "const.inc"

C INPUT/OUTPUT PARAMETERS

      integer istr, jstr, kstr, nip1, njp1, nkp1
      real    dtime, visc
      real    u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1),
     &        gam(nip1*njp1*nkp1), sc1(nip1*njp1*nkp1)
      real    rho(nip1*njp1*nkp1),Temperature(nip1*njp1*nkp1) 

      real    dip(nip1*njp1), djp(nip1*njp1), rp(nip1*njp1), dkw(nkp1)
      integer ijkmaxcfl(3), ijkmaxdfs(3)
      integer ijkmaxsc1(3), ijkminsc1(3)

      integer ijkmaxrho(3), ijkminrho(3) 

      integer ijkmaxT(3), ijkminT(3) 
      real    uvw_min(3), uvw_max(3)
      real    sc1_min, sc1_max

      real    rho_min, rho_max, T_min, T_max 
      real    cflmax, dfsmax

C LOCAL QUANTITIES

      integer i, j, k, joff, koff, jkoff, ij, ijk, jo, ko
      integer bigger, smaller
      integer imxcfl, jmxcfl, kmxcfl, imxdfs, jmxdfs, kmxdfs
      integer imaxsc1, jmaxsc1, kmaxsc1, iminsc1, jminsc1,kminsc1
      real    u_min, v_min, w_min,  u_max, v_max, w_max
      real    cfl, dfs
      real    dx, dy, dz, dx_inv, dy_inv, dz_inv

      integer imaxrho, jmaxrho, kmaxrho, iminrho, jminrho,kminrho 
      integer imaxT, jmaxT, kmaxT, iminT, jminT,kminT 

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C CALCULATE THE CFL NUMBER AT EVERY POINTS AND FIND THE MAXIMUM ONE

      cflmax = 0.0
      dfsmax = 0.0
      imxcfl = 0
      jmxcfl = 0
      kmxcfl = 0
      imxdfs = 0
      jmxdfs = 0
      kmxdfs = 0
      imaxsc1 = 0
      jmaxsc1 = 0
      kmaxsc1 = 0
      iminsc1 = 0
      jminsc1 = 0
      kminsc1 = 0

      imaxrho = 0 

      jmaxrho = 0

      kmaxrho = 0

      iminrho = 0

      jminrho = 0

      kminrho = 0 

c      imaxT = 0    

c      jmaxT = 0

c      kmaxT = 0

c      iminT = 0

c      jminT = 0

c      kminT = 0   
      u_max = -large_number
      v_max = -large_number
      w_max = -large_number
      sc1_max = -large_number

      rho_max = -large_number 

c      T_max = - large_number  
      u_min =  large_number
      v_min =  large_number
      w_min =  large_number
      sc1_min = large_number
      rho_min = large_number

c      T_min = large_number   



      jo = nip1
      ko = nip1*njp1

      do 10 k = kstr,nkp1-1
        koff = ko*(k-1)
        do 10 j = jstr,njp1-1
          joff  = jo*(j-1)
          jkoff = jo*(j-1) + koff
          do 10 i = istr,nip1-1
 
            ij  = i + joff
            ijk = i + jkoff
            dx = min(dip(ij-1 ),dip(ij))
            dy = min(djp(ij-jo),djp(ij))
            dz = rp(ij)*min(dkw(k),dkw(k-1))

            dx_inv = 1./dx
            dy_inv = 1./dy
            dz_inv = 1./dz
 
C MAXIMUM OF u, v, w

            u_max = max(u(ijk),u_max)
            v_max = max(v(ijk),v_max)
            w_max = max(w(ijk),w_max)
            if (sc1(ijk) .gt. sc1_max) then
               imaxsc1 = i
               jmaxsc1 = j
               kmaxsc1 = k
            endif
            sc1_max = max(sc1(ijk),sc1_max)

C MAXIMUM OF rho & T
            if (rho(ijk) .gt. rho_max) then
               imaxrho = i
               jmaxrho = j
               kmaxrho = k
            endif
            rho_max = max(rho(ijk),rho_max)

c            if (Temperature(ijk) .gt. T_max) then
c               imaxT = i
c               jmaxT = j
c               kmaxT = k
c            endif
c            T_max = max(Temperature(ijk),T_max)


C MINIMUM OF u, v, w

            u_min = min(u(ijk),u_min)
            v_min = min(v(ijk),v_min)
            w_min = min(w(ijk),w_min)
            if (sc1(ijk) .lt. sc1_min) then
               iminsc1 = i
               jminsc1 = j
               kminsc1 = k
            endif
            sc1_min = min(sc1(ijk),sc1_min)

C MINIMUM OF rho & T
            if (rho(ijk) .lt. rho_min) then
               iminrho = i
               jminrho = j
               kminrho = k
            endif
            rho_min = min(rho(ijk),rho_min)

c            if (Temperature(ijk) .lt. T_min) then
c               iminT = i
c               jminT = j
c               kminT = k
c            endif
c            T_min = min(Temperature(ijk),T_min)


C CFL NUMBER

            cfl = dtime*
     &        (abs(u(ijk))*dx_inv + abs(v(ijk))*dy_inv +
     &         abs(w(ijk))*dz_inv)

            bigger  = int((sign(1.0,cflmax-cfl) + 1.0)*0.5 + 0.5)
            smaller = 1 - bigger

            imxcfl = imxcfl*bigger + i*smaller
            jmxcfl = jmxcfl*bigger + j*smaller
            kmxcfl = kmxcfl*bigger + k*smaller
            cflmax = max(cflmax,cfl)

C THE NUMBER FOR JUDGING DIFFUSIVE STABILITY

cwk            dfs = dtime*(gam(ijk)+visc)*
cwk     &        (dx_inv*dx_inv + dy_inv*dy_inv + dz_inv*dz_inv)
            dfs = dtime*(gam(ijk)+visc/rho(ijk))*
     &        (dx_inv*dx_inv + dy_inv*dy_inv + dz_inv*dz_inv)  


            bigger  = int((sign(1.0,dfsmax-dfs) + 1.0)*0.5 + 0.5)
            smaller = 1 - bigger

            imxdfs = imxdfs*bigger + i*smaller
            jmxdfs = jmxdfs*bigger + j*smaller
            kmxdfs = kmxdfs*bigger + k*smaller
            dfsmax = max(dfsmax,dfs)

 20       continue
 15     continue
 10   continue

      ijkmaxcfl(1) = imxcfl
      ijkmaxcfl(2) = jmxcfl
      ijkmaxcfl(3) = kmxcfl
      ijkmaxdfs(1) = imxdfs
      ijkmaxdfs(2) = jmxdfs
      ijkmaxdfs(3) = kmxdfs
      ijkmaxsc1(1) = imaxsc1
      ijkmaxsc1(2) = jmaxsc1
      ijkmaxsc1(3) = kmaxsc1
      ijkminsc1(1) = iminsc1
      ijkminsc1(2) = jminsc1
      ijkminsc1(3) = kminsc1

      uvw_min(1) = u_min
      uvw_min(2) = v_min
      uvw_min(3) = w_min
      uvw_max(1) = u_max
      uvw_max(2) = v_max
      uvw_max(3) = w_max


      ijkmaxrho(1) = imaxrho
      ijkmaxrho(2) = jmaxrho
      ijkmaxrho(3) = kmaxrho
      ijkminrho(1) = iminrho
      ijkminrho(2) = jminrho
      ijkminrho(3) = kminrho

c      ijkmaxT(1) = imaxT
c      ijkmaxT(2) = jmaxT
c      ijkmaxT(3) = kmaxT
c      ijkminT(1) = iminT
c      ijkminT(2) = jminT
c      ijkminT(3) = kminT


      return
      end
