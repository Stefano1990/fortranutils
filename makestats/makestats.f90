program statistics

implicit none

  integer :: i,j,k,maxbl,n,m,jlow,nib,njb
  parameter (maxbl = 128)
  integer :: ni,nj,t,tmax,nk,off1,off2,row,bpr
  parameter (ni = 1024, nj=256,nk=134,nib=32,njb=64,bpr=32)
  double precision :: x(1:nib,1:njb,1:maxbl), y(1:nib,1:njb,1:maxbl)
  double precision :: u(1:nib,1:njb,1:maxbl), v(1:nib,1:njb,1:maxbl),w(1:nib,1:njb,1:maxbl)
  double precision:: sc1(1:nib,1:njb,1:maxbl),p(1:nib,1:njb,1:maxbl),rho(1:nib,1:njb,1:maxbl)
  double precision :: sc1var(1:nib,1:njb,1:maxbl),evisc(1:nib,1:njb,1:maxbl)
  double precision :: wu(1:nib,1:njb,1:maxbl),vw(1:nib,1:njb,1:maxbl)     
  double precision :: xp(1:ni,1:nj),yp(1:ni,1:nj),up(1:ni,1:nj),vp(1:ni,1:nj)
  double precision :: wp(1:ni,1:nj),sc1p(1:ni,1:nj),pp(1:ni,1:nj)
  double precision :: rhop(1:ni,1:nj),sc1varp(1:ni,1:nj)
  double precision :: uu(1:nib,1:njb,1:maxbl), vv(1:nib,1:njb,1:maxbl),ww(1:nib,1:njb,1:maxbl)
  double precision :: uv(1:nib,1:njb,1:maxbl),ff(1:nib,1:njb,1:maxbl),rr(1:nib,1:njb,1:maxbl)
  double precision :: uup(1:ni,1:nj),vvp(1:ni,1:nj)
  double precision :: wwp(1:ni,1:nj),uvp(1:ni,1:nj),ffp(1:ni,1:nj)
  double precision :: rrp(1:ni,1:nj)
  double precision :: uufluc(1:ni,1:nj),vvfluc(1:ni,1:nj)
  double precision :: wwfluc(1:ni,1:nj),uvfluc(1:ni,1:nj),fffluc(1:ni,1:nj)
  double precision :: rrfluc(1:ni,1:nj),temp, dudy
  double precision :: dp(1:ni), du(1:ni)
  double precision :: Rl(1:ni), dw(1:ni),yu(1:ni,1:nj)
  
!open(unit=14,file='mean.dat')
!open(unit=15,file='stress.dat')
open(unit=16,file='statistics.dat')

do i=1, maxbl
 open(unit=5000+i)
 open(unit=6000+i)
 read(5000+i,*)
 read(5000+i,*)
 read(6000+i,*)
 read(6000+i,*)
enddo

!read(14,*)
!read(14,*)
!read(15,*)
!read(15,*)

! Does all the rows in one go.

do n=1,maxbl
 read(5000+n,*)
 write(*,*) 'reading 5000 blocks for bl nr:',n
 do j=1,njb
  do i=1,nib
   read(5000+n,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
      sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
  enddo
 enddo
 
 read(6000+n,*)
 write(*,*) 'reading 6000 blocks for bl nr:',n
 do j=1,njb
  do i=1,nib
   read(6000+n,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n), &
      wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
  enddo
 enddo
enddo


do row=1,4
do n=bpr*row-(bpr-1),row*bpr
 do j=1,njb
  do i=1,nib
   xp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = x(i,j,n)
   yp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = y(i,j,n)
   up(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = u(i,j,n)/rho(i,j,n)
   vp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = v(i,j,n)/rho(i,j,n)
   wp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = w(i,j,n)/rho(i,j,n)
   sc1p(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = sc1(i,j,n)
   pp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = p(i,j,n)
   rhop(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = rho(i,j,n)
   uup(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = uu(i,j,n)/rho(i,j,n)
   vvp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = vv(i,j,n)/rho(i,j,n)
   wwp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = ww(i,j,n)/rho(i,j,n)
   uvp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = uv(i,j,n)/rho(i,j,n)
   ffp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = ff(i,j,n)
   rrp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = rr(i,j,n)
   sc1varp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb) = sc1var(i,j,n)

   if(i+(n-(bpr*row-(bpr-1)))*nib == 1024 .and. j+(row-1)*njb == 128) then ! i=j=1
    write(*,*) 'i=1,j=256;x=',xp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb),'y=',yp(i+(n-(bpr*row-(bpr-1)))*nib,j+(row-1)*njb)
   endif
  enddo
 enddo
 enddo
enddo
!do n=17,32
! do j=1,njb
!  do i=1,nib
!   xp(i+(n-17)*ni,j+njb) = x(i,j,n)
!   yp(i+(n-17)*ni,j+njb) = y(i,j,n)
!   up(i+(n-17)*ni,j+njb) = u(i,j,n)/rho(i,j,n)
!   vp(i+(n-17)*ni,j+njb) = v(i,j,n)/rho(i,j,n)
!   wp(i+(n-17)*ni,j+njb) = w(i,j,n)/rho(i,j,n)
!   sc1p(i+(n-17)*ni,j+njb) = sc1(i,j,n)
!   pp(i+(n-17)*ni,j+njb) = p(i,j,n)
!   rhop(i+(n-17)*ni,j+njb) = rho(i,j,n)
!   uup(i+(n-17)*ni,j+njb) = uu(i,j,n)/rho(i,j,n)
!   vvp(i+(n-17)*ni,j+njb) = vv(i,j,n)/rho(i,j,n)
!   wwp(i+(n-17)*ni,j+njb) = ww(i,j,n)/rho(i,j,n)
!   uvp(i+(n-17)*ni,j+njb) = uv(i,j,n)/rho(i,j,n)
!   ffp(i+(n-17)*ni,j+njb) = ff(i,j,n)
!   rrp(i+(n-17)*ni,j+njb) = rr(i,j,n)
!   sc1varp(i+(n-17)*ni,j+njb) = sc1var(i,j,n)
!  enddo
! enddo
!enddo

!do n=1,16
! read(5000+n,*)
! write(*,*) 'reading 5000 blocks for bl nr:',n
! do j=1,njb
!  do i=1,nib
!   read(5000+n,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
!      sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=1,16
! read(6000+n,*)
! write(*,*) 'reading 6000 blocks for bl nr:',n
! do j=1,njb
!  do i=1,nib
!   read(6000+n,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n), &
!      wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=1,16
! do j=1,njb
!  do i=1,nib
!   xp(i+(n-1)*ni,j) = x(i,j,n)
!   yp(i+(n-1)*ni,j) = y(i,j,n)
!   up(i+(n-1)*ni,j) = u(i,j,n)/rho(i,j,n)
!   vp(i+(n-1)*ni,j) = v(i,j,n)/rho(i,j,n)
!   wp(i+(n-1)*ni,j) = w(i,j,n)/rho(i,j,n)
!   sc1p(i+(n-1)*ni,j) = sc1(i,j,n)
!   pp(i+(n-1)*ni,j) = p(i,j,n)
!   rhop(i+(n-1)*ni,j) = rho(i,j,n)
!   uup(i+(n-1)*ni,j) = uu(i,j,n)/rho(i,j,n)
!   vvp(i+(n-1)*ni,j) = vv(i,j,n)/rho(i,j,n)
!   wwp(i+(n-1)*ni,j) = ww(i,j,n)/rho(i,j,n)
!   uvp(i+(n-1)*ni,j) = uv(i,j,n)/rho(i,j,n)
!   ffp(i+(n-1)*ni,j) = ff(i,j,n)
!   rrp(i+(n-1)*ni,j) = rr(i,j,n)
!   sc1varp(i+(n-1)*ni,j) = sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!! Second row
!do n=17,32
! read(5000+n,*)
! do j=1,njb
!  do i=1,nib
!   read(5000+n,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
!      sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=17,32
! read(6000+n,*)
! do j=1,njb
!  do i=1,nib
!   read(6000+n,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n), &
!      wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=17,32
! do j=1,njb
!  do i=1,nib
!   xp(i+(n-17)*ni,j+njb) = x(i,j,n)
!   yp(i+(n-17)*ni,j+njb) = y(i,j,n)
!   up(i+(n-17)*ni,j+njb) = u(i,j,n)/rho(i,j,n)
!   vp(i+(n-17)*ni,j+njb) = v(i,j,n)/rho(i,j,n)
!   wp(i+(n-17)*ni,j+njb) = w(i,j,n)/rho(i,j,n)
!   sc1p(i+(n-17)*ni,j+njb) = sc1(i,j,n)
!   pp(i+(n-17)*ni,j+njb) = p(i,j,n)
!   rhop(i+(n-17)*ni,j+njb) = rho(i,j,n)
!   uup(i+(n-17)*ni,j+njb) = uu(i,j,n)/rho(i,j,n)
!   vvp(i+(n-17)*ni,j+njb) = vv(i,j,n)/rho(i,j,n)
!   wwp(i+(n-17)*ni,j+njb) = ww(i,j,n)/rho(i,j,n)
!   uvp(i+(n-17)*ni,j+njb) = uv(i,j,n)/rho(i,j,n)
!   ffp(i+(n-17)*ni,j+njb) = ff(i,j,n)
!   rrp(i+(n-17)*ni,j+njb) = rr(i,j,n)
!   sc1varp(i+(n-17)*ni,j+njb) = sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!! third row
!do n=33,48
! read(5000+n,*)
! do j=1,njb
!  do i=1,nib
!   read(5000+n,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
!      sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=33,48
! read(6000+n,*)
! do j=1,njb
!  do i=1,nib
!   read(6000+n,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n), &
!      wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=33,48
! do j=1,njb
!  do i=1,nib
!   xp(i+(n-33)*ni,j+2*njb) = x(i,j,n)
!   yp(i+(n-33)*ni,j+2*njb) = y(i,j,n)
!   up(i+(n-33)*ni,j+2*njb) = u(i,j,n)/rho(i,j,n)
!   vp(i+(n-33)*ni,j+2*njb) = v(i,j,n)/rho(i,j,n)
!   wp(i+(n-33)*ni,j+2*njb) = w(i,j,n)/rho(i,j,n)
!   sc1p(i+(n-33)*ni,j+2*njb) = sc1(i,j,n)
!   pp(i+(n-33)*ni,j+2*njb) = p(i,j,n)
!   rhop(i+(n-33)*ni,j+2*njb) = rho(i,j,n)
!   uup(i+(n-33)*ni,j+2*njb) = uu(i,j,n)/rho(i,j,n)
!   vvp(i+(n-33)*ni,j+2*njb) = vv(i,j,n)/rho(i,j,n)
!   wwp(i+(n-33)*ni,j+2*njb) = ww(i,j,n)/rho(i,j,n)
!   uvp(i+(n-33)*ni,j+2*njb) = uv(i,j,n)/rho(i,j,n)
!   ffp(i+(n-33)*ni,j+2*njb) = ff(i,j,n)
!   rrp(i+(n-33)*ni,j+2*njb) = rr(i,j,n)
!   sc1varp(i+(n-33)*ni,j+2*njb) = sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!! Fourth row
!do n=49,64
! read(5000+n,*)
! do j=1,njb
!  do i=1,nib
!   read(5000+n,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
!      sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=49,64
! read(6000+n,*)
! do j=1,njb
!  do i=1,nib
!   read(6000+n,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n), &
!      wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!do n=49,64
! do j=1,njb
!  do i=1,nib
!   xp(i+(n-49)*ni,j+3*njb) = x(i,j,n)
!   yp(i+(n-49)*ni,j+3*njb) = y(i,j,n)
!   up(i+(n-49)*ni,j+3*njb) = u(i,j,n)/rho(i,j,n)
!   vp(i+(n-49)*ni,j+3*njb) = v(i,j,n)/rho(i,j,n)
!   wp(i+(n-49)*ni,j+3*njb) = w(i,j,n)/rho(i,j,n)
!   sc1p(i+(n-49)*ni,j+3*njb) = sc1(i,j,n)
!   pp(i+(n-49)*ni,j+3*njb) = p(i,j,n)
!   rhop(i+(n-49)*ni,j+3*njb) = rho(i,j,n)
!   uup(i+(n-49)*ni,j+3*njb) = uu(i,j,n)/rho(i,j,n)
!   vvp(i+(n-49)*ni,j+3*njb) = vv(i,j,n)/rho(i,j,n)
!   wwp(i+(n-49)*ni,j+3*njb) = ww(i,j,n)/rho(i,j,n)
!   uvp(i+(n-49)*ni,j+3*njb) = uv(i,j,n)/rho(i,j,n)
!   ffp(i+(n-49)*ni,j+3*njb) = ff(i,j,n)
!   rrp(i+(n-49)*ni,j+3*njb) = rr(i,j,n)
!   sc1varp(i+(n-49)*ni,j+3*njb) = sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!do j=1,njb
! do i=1,nib
!  do n=1,maxbl
!   if(n >= 1 .and. n <= 16) then
!    off1 = 1
!    off2 = 0
!   endif
!   if(n >= 17 .and. n <= 32) then
!    off1 = 17
!    off2 = njb
!   endif
!   if(n >= 33 .and. n <= 48) then
!    off1 = 33
!    off2 = 2*njb
!   endif
!   if(n >= 49 .and. n <= 64) then
!    off1 = 49
!    off2 = 3*njb
!   endif
!   xp(i+(n-off1)*nib,j+off2) = x(i,j,n)
!   yp(i+(n-off1)*nib,j+off2) = y(i,j,n)
!   up(i+(n-off1)*nib,j+off2) = u(i,j,n)/rho(i,j,n)
!   vp(i+(n-off1)*nib,j+off2) = v(i,j,n)/rho(i,j,n)
!   wp(i+(n-off1)*nib,j+off2) = w(i,j,n)/rho(i,j,n)
!   sc1p(i+(n-off1)*nib,j+off2) = sc1(i,j,n)
!   pp(i+(n-off1)*nib,j+off2) = p(i,j,n)
!   rhop(i+(n-off1)*nib,j+off2) = rho(i,j,n)
!   uup(i+(n-off1)*nib,j+off2) = uu(i,j,n)/rho(i,j,n)
!   vvp(i+(n-off1)*nib,j+off2) = vv(i,j,n)/rho(i,j,n)
!   wwp(i+(n-off1)*nib,j+off2) = ww(i,j,n)/rho(i,j,n)
!   uvp(i+(n-off1)*nib,j+off2) = uv(i,j,n)/rho(i,j,n)
!   ffp(i+(n-off1)*nib,j+off2) = ff(i,j,n)
!   rrp(i+(n-off1)*nib,j+off2) = rr(i,j,n)
!   sc1varp(i+(n-off1)*nib,j+off2) = sc1var(i,j,n)
!   if(vp(i+(n-off1)*nib,j+off2) /= vp(i+(n-off1)*ni,j+off2)) then ! NaN
!    write(*,*) 'NaN at:',i+(n-off1)*nib,j+off2,n
!   endif
!  enddo
! enddo
!enddo
!
! do k=1,maxbl 
!   m=5000+k
!   n=k
!   read(m,*)
!
!   do j=1,nj
!    do i=1, ni
!     read(m,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
!              sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
!   enddo
!  enddo
! enddo
!
!do k=1, maxbl
! m=6000+k
! n=k
! read(m,*)
! do j=1,nj
!  do i=1,ni
!  read(m,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n),wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
!  enddo
! enddo
!enddo
!
!
!do j=1,nj
! do i=1,ni 
! do n=1,maxbl
!  xp(i+(n-1)*ni) = x(i,n)
!  up(i+(n-1)*ni,j) = u(i,j,n)/rho(i,j,n)
!  vp(i+(n-1)*ni,j) = v(i,j,n)/rho(i,j,n)
!  wp(i+(n-1)*ni,j) = w(i,j,n)/rho(i,j,n)
!  sc1p(i+(n-1)*ni,j) = sc1(i,j,n)
!  pp(i+(n-1)*ni,j) = p(i,j,n)
!  rhop(i+(n-1)*ni,j) = rho(i,j,n)
!  uup(i+(n-1)*ni,j) = uu(i,j,n)/rho(i,j,n)
!  vvp(i+(n-1)*ni,j) = vv(i,j,n)/rho(i,j,n)
!  wwp(i+(n-1)*ni,j) = ww(i,j,n)/rho(i,j,n)
!  uvp(i+(n-1)*ni,j) = uv(i,j,n)/rho(i,j,n)
!  ffp(i+(n-1)*ni,j) = ff(i,j,n)
!  rrp(i+(n-1)*ni,j) = rr(i,j,n)
!  sc1varp(i+(n-1)*ni,j) = sc1var(i,j,n)
! enddo
! enddo
!enddo
!
do i=1,ni
 Rl(i) = (up(i,nj-3) - up(i,2)) / (up(i,nj-3) + up(i,2))
 du(i) = up(i,nj-3) - up(i,2)
 dp(i) = pp(i,nj-3) - pp(i,2)
enddo

do j=1,nj
 do i=1,ni
  uufluc(i,j) = uup(i,j) - up(i,j)**2
  vvfluc(i,j) = vvp(i,j) - vp(i,j)**2
  wwfluc(i,j) = wwp(i,j) - wp(i,j)**2
  uvfluc(i,j) = (uup(i,j) - up(i,j)**2)*(vvp(i,j) - vp(i,j)**2)
  fffluc(i,j) = ffp(i,j) - sc1p(i,j)**2
  rrfluc(i,j) = rrp(i,j) - rhop(i,j)**2 
 enddo
enddo

write(16,*) 'VARIABLES = "x" "y" "u" "v" "w" "f" "p" "rho" "uu" "vv" "ww" "uv" "ff" "rr" "scalvar"'
write(16,*) 'ZONE I=', ni,', J=', nj
do j=1,nj
 do i=1,ni
  write(16,*) xp(i,j), yp(i,j), up(i,j),vp(i,j),wp(i,j),sc1p(i,j), pp(i,j), rhop(i,j), &
      uufluc(i,j), vvfluc(i,j), wwfluc(i,j), uvfluc(i,j), fffluc(i,j), rrfluc(i,j),sc1varp(i,j)
            
 enddo
enddo
do i=1,ni
do j=1, nj-1
 yu(i,j) = 0.5*(yp(i,j) + yp(i,j+1))
!write(*,*) j, yu(j)
enddo
enddo

do i = 1, ni
temp = 0.0
do j=  2, nj-1
dudy = (up(i,j) - up(i,j-1)) / (yu(i,j) - yu(i,j-1))
if(dudy .gt. temp) then
temp = dudy
jlow = j
endif
enddo
if(i.ge.189.and. i.le.193) write(*,*) i, temp,jlow
dw(i) = (up(i,nj-2) - up(i,2)) / temp
enddo
!
!
open(unit=17,file='dp.dat')
write(17,*) 'VARIABLES = "x" "Rl" "du" "dp" "vthick"'
do i=1, ni
write(17,*) xp(i,1), Rl(i), du(i),dp(i),dw(i)
enddo
!

end program
