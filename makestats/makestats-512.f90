program statistics

implicit none

  integer :: i,j,k,maxbl,n,m,jlow,nib,njb
  parameter (maxbl = 64)
  integer :: ni,nj,t,tmax,nk,off1,off2,row,bpr,rows
  integer :: ig,jg
  parameter (ni = 512, nj=128,nk=128,nib=32,njb=32,bpr=16,rows=4)
  double precision :: x(1:nib,1:njb,1:maxbl), y(1:nib,1:njb,1:maxbl)
  double precision :: u(1:nib,1:njb,1:maxbl), v(1:nib,1:njb,1:maxbl),w(1:nib,1:njb,1:maxbl)
  double precision:: sc1(1:nib,1:njb,1:maxbl),p(1:nib,1:njb,1:maxbl),rho(1:nib,1:njb,1:maxbl)
  double precision :: sc1var(1:nib,1:njb,1:maxbl)
  double precision :: wu(1:nib,1:njb,1:maxbl),vw(1:nib,1:njb,1:maxbl)     
  double precision :: xp(1:ni,1:nj),yp(1:ni,1:nj),up(1:ni,1:nj),vp(1:ni,1:nj)
  double precision :: wp(1:ni,1:nj),sc1p(1:ni,1:nj),pp(1:ni,1:nj)
  double precision :: rhop(1:ni,1:nj),sc1varp(1:ni,1:nj)
  double precision :: uu(1:nib,1:njb,1:maxbl), vv(1:nib,1:njb,1:maxbl),ww(1:nib,1:njb,1:maxbl)
  double precision :: uv(1:nib,1:njb,1:maxbl),ff(1:nib,1:njb,1:maxbl),rr(1:nib,1:njb,1:maxbl)
  double precision :: uup(1:ni,1:nj),vvp(1:ni,1:nj)
  double precision :: wwp(1:ni,1:nj),uvp(1:ni,1:nj),ffp(1:ni,1:nj)
  double precision :: rrp(1:ni,1:nj),evisc(1:nib,1:njb,1:maxbl)
  double precision :: uufluc(1:ni,1:nj),vvfluc(1:ni,1:nj)
  double precision :: wwfluc(1:ni,1:nj),uvfluc(1:ni,1:nj),fffluc(1:ni,1:nj)
  double precision :: rrfluc(1:ni,1:nj),temp, dudy
  double precision :: dp(1:ni), du(1:ni),eviscp(1:ni,1:nj)
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

! first row
do n=1,maxbl
read(6000+n,*)
! write(*,*) 'reading 6000 blocks for bl nr:',n
 do j=1,njb
  do i=1,nib
   read(6000+n,*) x(i,j,n), y(i,j,n), uu(i,j,n),vv(i,j,n),ww(i,j,n),uv(i,j,n), &
      wu(i,j,n),vw(i,j,n),ff(i,j,n),rr(i,j,n),sc1var(i,j,n)
  enddo
 enddo

read(5000+n,*)
! write(*,*) 'reading 5000 blocks for bl nr:',n
 do j=1,njb
  do i=1,nib
   read(5000+n,*) x(i,j,n), y(i,j,n), u(i,j,n),v(i,j,n),w(i,j,n),&
      sc1(i,j,n),p(i,j,n),rho(i,j,n),evisc(i,j,n)
  enddo
 enddo
enddo

do n=1,16
 do j=1,njb
  do i=1,nib
   ig = i+(n-1)*nib
   jg = j
   xp(ig,jg) = x(i,j,n)
   yp(ig,jg) = y(i,j,n)
   up(ig,jg) = u(i,j,n)/rho(i,j,n)
   vp(ig,jg) = v(i,j,n)/rho(i,j,n)
   wp(ig,jg) = w(i,j,n)/rho(i,j,n)
   sc1p(ig,jg) = sc1(i,j,n)
   pp(ig,jg) = p(i,j,n)
   rhop(ig,jg) = rho(i,j,n)
   uup(ig,jg) = uu(i,j,n)/rho(i,j,n)
   vvp(ig,jg) = vv(i,j,n)/rho(i,j,n)
   wwp(ig,jg) = ww(i,j,n)/rho(i,j,n)
   uvp(ig,jg) = uv(i,j,n)/rho(i,j,n)
   ffp(ig,jg) = ff(i,j,n)
   rrp(ig,jg) = rr(i,j,n)
   sc1varp(ig,jg) = sc1var(i,j,n)
   eviscp(ig,jg) = evisc(i,j,n)
  enddo
 enddo
enddo

do n=17,32
 do j=1,njb
  do i=1,nib
   ig = i+(n-17)*nib
   jg = j+njb
   xp(ig,jg) = x(i,j,n)
   yp(ig,jg) = y(i,j,n)
   up(ig,jg) = u(i,j,n)/rho(i,j,n)
   vp(ig,jg) = v(i,j,n)/rho(i,j,n)
   wp(ig,jg) = w(i,j,n)/rho(i,j,n)
   sc1p(ig,jg) = sc1(i,j,n)
   pp(ig,jg) = p(i,j,n)
   rhop(ig,jg) = rho(i,j,n)
   uup(ig,jg) = uu(i,j,n)/rho(i,j,n)
   vvp(ig,jg) = vv(i,j,n)/rho(i,j,n)
   wwp(ig,jg) = ww(i,j,n)/rho(i,j,n)
   uvp(ig,jg) = uv(i,j,n)/rho(i,j,n)
   ffp(ig,jg) = ff(i,j,n)
   rrp(ig,jg) = rr(i,j,n)
   sc1varp(ig,jg) = sc1var(i,j,n)
   eviscp(ig,jg) = evisc(i,j,n)
  enddo
 enddo
enddo

do n=33,48
 do j=1,njb
  do i=1,nib
   ig = i+(n-33)*nib
   jg = j+njb*2
   xp(ig,jg) = x(i,j,n)
   yp(ig,jg) = y(i,j,n)
   up(ig,jg) = u(i,j,n)/rho(i,j,n)
   vp(ig,jg) = v(i,j,n)/rho(i,j,n)
   wp(ig,jg) = w(i,j,n)/rho(i,j,n)
   sc1p(ig,jg) = sc1(i,j,n)
   pp(ig,jg) = p(i,j,n)
   rhop(ig,jg) = rho(i,j,n)
   uup(ig,jg) = uu(i,j,n)/rho(i,j,n)
   vvp(ig,jg) = vv(i,j,n)/rho(i,j,n)
   wwp(ig,jg) = ww(i,j,n)/rho(i,j,n)
   uvp(ig,jg) = uv(i,j,n)/rho(i,j,n)
   ffp(ig,jg) = ff(i,j,n)
   rrp(ig,jg) = rr(i,j,n)
   sc1varp(ig,jg) = sc1var(i,j,n)
   eviscp(ig,jg) = evisc(i,j,n)
  enddo
 enddo
enddo

do n=49,64
 do j=1,njb
  do i=1,nib
   ig = i+(n-49)*nib
   jg = j+njb*3
   xp(ig,jg) = x(i,j,n)
   yp(ig,jg) = y(i,j,n)
   up(ig,jg) = u(i,j,n)/rho(i,j,n)
   vp(ig,jg) = v(i,j,n)/rho(i,j,n)
   wp(ig,jg) = w(i,j,n)/rho(i,j,n)
   sc1p(ig,jg) = sc1(i,j,n)
   pp(ig,jg) = p(i,j,n)
   rhop(ig,jg) = rho(i,j,n)
   uup(ig,jg) = uu(i,j,n)/rho(i,j,n)
   vvp(ig,jg) = vv(i,j,n)/rho(i,j,n)
   wwp(ig,jg) = ww(i,j,n)/rho(i,j,n)
   uvp(ig,jg) = uv(i,j,n)/rho(i,j,n)
   ffp(ig,jg) = ff(i,j,n)
   rrp(ig,jg) = rr(i,j,n)
   sc1varp(ig,jg) = sc1var(i,j,n)
   eviscp(ig,jg) = evisc(i,j,n)
  enddo
 enddo
enddo


write(*,*) xp(33,1),yp(33,1)

do i=1,1
 do j=1,1
  do n=33,33
 !  write(*,*) i,j,xp(i,j),yp(i,j)
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

write(16,*) 'VARIABLES = "x" "y" "u" "v" "w" "f" "p" "rho" "uu" "vv" "ww" "uv" "ff" "rr" "sc1var" "evisc"'
write(16,*) 'ZONE I=', ni,', J=', nj
do j=1,nj
 do i=1,ni
  write(16,*) xp(i,j), yp(i,j), up(i,j),vp(i,j),wp(i,j),sc1p(i,j), pp(i,j), rhop(i,j), &
      uufluc(i,j), vvfluc(i,j), wwfluc(i,j), uvfluc(i,j), fffluc(i,j), rrfluc(i,j), &
      sc1varp(i,j),eviscp(i,j)
            
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
!if(i.ge.189.and. i.le.193) write(*,*) i, temp,jlow
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
