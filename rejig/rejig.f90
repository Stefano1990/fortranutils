program flowrearrange


implicit none

  integer :: i,j,k,maxbl,n,m,recycbls,row,rows,bpr
  parameter (maxbl = 128, recycbls = 16,bpr=32,rows=4)
  integer :: ni,nj,nj1,nj2,t,tmax,nk,tmin
  parameter (ni = 32, nj=256,nj1=64,nj2=100,nk=128)
   double precision :: x(1:ni+1,1:nj+1,1:maxbl), y(1:ni+1,1:nj+1,1:maxbl)
  double precision :: u(1:ni+1,1:nj+1,1:maxbl), v(1:ni+1,1:nj+1,1:maxbl),w(1:ni+1,1:nj+1,1:maxbl)
  double precision :: sc1(1:ni+1,1:nj+1,1:maxbl),drho(1:ni+1,1:nj+1,1:maxbl)
  double precision :: xp(1:maxbl*ni,1:nj), up(1:maxbl*ni,1:nj),vp(1:maxbl*ni,1:nj)
  double precision :: yp(1:maxbl*ni,1:nj),d2rho(1:ni+1,1:nj+1,1:maxbl)
  double precision :: wp(1:maxbl*ni,1:nj),sc1p(1:maxbl*ni,1:nj),drhop(1:maxbl*ni,1:nj)
  double precision :: d2rhop(1:maxbl*ni,1:nj),d2rhopy(1:maxbl*ni,1:nk)
  double precision :: xdum(1:ni,1:maxbl), z(1:nk,1:maxbl),drhopy(1:maxbl*ni,1:nk)
  double precision :: drhoy(1:ni,1:nk,1:maxbl),d2rhoy(1:ni,1:nk,1:maxbl), uc(1:ni*maxbl) 
  
open(unit=16, file='flowviz.dat')
write(16,*) 'VARIABLES = "x" "y" "u" "u-uc" "v" "w" "sc1" '

  print *, 'tmin = ?'
  read *, tmin
  
  print *, 'tmax = ?'
  read *, tmax
   
   
!print *, 'Type in high- and low-speed stream densities'
!read *, rhohi, rholow



!  rhoparam = (sqrt(rhohi) - sqrt(rholow)) / (sqrt(rhohi) + sqrt(rholow))

do i=1, maxbl
open(unit=3000+i)
read(3000+i,*)
enddo

t=1
do while (t .le. tmax)
write(*,*) t


!get the first row of blocks into their part of the zone
do n=1, maxbl
   m = 3000+n
read(m,*)
  do j=1,nj1+1
  do i=1, ni+1
   read(m,*) x(i,j,n), y(i,j,n),u(i,j,n),v(i,j,n),w(i,j,n),sc1(i,j,n),drho(i,j,n),d2rho(i,j,n)
   enddo
  enddo
 enddo

do row=1,rows
do n=bpr*row-(bpr-1),maxbl
 do j=1,nj1
  do i=1,ni
   xp(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = x(i,j,n)
   yp(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = y(i,j,n)
   up(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = u(i,j,n)
   vp(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = v(i,j,n)
   wp(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = w(i,j,n)
   sc1p(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = sc1(i,j,n)
   drhop(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = drho(i,j,n)
   d2rhop(i+(n-(bpr*row-(bpr-1)))*ni,j+(row-1)*nj1) = d2rho(i,j,n)
  enddo
 enddo
 enddo
enddo


!do i=1,16*ni
! do j = nj1+1, 2*nj
!  yp(i,j) = y(i,j-nj1,17)
! enddo
! do j=2*nj1+1,3*nj
!  yp(i,j) = y(i,j-2*nj1,33)
! enddo
! do j=3*nj1+1, 4*nj
!  yp(i,j) = y(i,j-3*nj1,49)
! enddo
!enddo
if(t >= tmin .and. t <= tmax) then
do i=1, bpr*ni
 uc(i) = 0.5*(up(i,nj-1) + up(i,2))
enddo

write(16,*) 'ZONE I=', ni*bpr,', J=', nj

do j=1, nj
do i=1, bpr*ni
!  lambda(i) = (up(i,nj-1) - up(i,2))/(up(i,nj-2) + up(i,3))
!  uc(i) = 0.5*(up(i,2) + up(1,nj-1))*(1+rhoparam*lambda(i))
write(16,*) xp(i,j), yp(i,j), up(i,j), up(i,j) - uc(i),vp(i,j), wp(i,j), sc1p(i,j)
enddo
enddo
endif

t=t+1
enddo

end
