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
            
            dudi = (u(ijke)-u(ijk))/di + vp*hij(i,j)
            dvdi = (vef    -vwf   )/di - up*hij(i,j)
            
            dudj = (unf    -usf   )/dj - vp*hji(i,j)
            dvdj = (v(ijkn)-v(ijk))/dj + up*hji(i,j)

C THE REST WITH RADIUS
C  - DIFFERENCE W/R => NO STRESS FOR SOLID BODY ROTATION ???????????????

            if (gotrad) then
            
              radp = rp(i  ,j  )
              radw = rp(i-1,j  )
              rade = rp(i+1,j  )
              rads = rp(i  ,j-1)
              radn = rp(i  ,j+1)
              
              dk = radp*dkp(i,j)
              
              hki  = (ru(i+1,j  )-ru(i,j))/ (radp*di)
              hkj  = (rv(i  ,j+1)-rv(i,j))/ (radp*dj)

              wdrw = 0.5*(w(ijkw) + w(ijkwr))/radw
              wdre = 0.5*(w(ijke) + w(ijker))/rade
              wdrs = 0.5*(w(ijks) + w(ijksr))/rads
              wdrn = 0.5*(w(ijkn) + w(ijknr))/radn
              
              dwdi = radp*(wdre - wdrw)/(diu(i,j)+diu(i+1,j  ))
              dwdj = radp*(wdrn - wdrs)/(djv(i,j)+djv(i  ,j+1))

              dudk = (urf    -ulf   )/dk
              dvdk = (vrf    -vlf   )/dk
              dwdk = (w(ijke)-w(ijk))/dk + up*hki + vp*hkj

C THE REST WITHOUT RADIUS

            else
              
              dk = dkp(i,j)

              wwf = 0.25*(w(ijk) + w(ijkr) + w(ijkw) + w(ijkwr))
              wef = 0.25*(w(ijk) + w(ijkr) + w(ijke) + w(ijker))
              wsf = 0.25*(w(ijk) + w(ijkr) + w(ijks) + w(ijksr))
              wnf = 0.25*(w(ijk) + w(ijkr) + w(ijkn) + w(ijknr))
              
              dwdi = (wef    -wwf   )/di
              dwdj = (wnf    -wsf   )/dj
              
              dudk = (urf    -ulf   )/dk
              dvdk = (vrf    -vlf   )/dk
              dwdk = (w(ijke)-w(ijk))/dk
              
            end if
