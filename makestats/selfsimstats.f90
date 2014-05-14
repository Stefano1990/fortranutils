program selfstats

implicit none

integer :: i,j,k,ni,nj,block1ni,block2ni,nj1,nj2,jhi,jlo
parameter(ni = 1024, nj=240,nj1=120,nj2=100,block1ni=768,block2ni=256)
double precision :: x(1:ni),y(1:nj)
double precision :: u(1:ni,1:nj),v(1:ni,1:nj),w(1:ni,1:nj)
double precision :: xi(1:ni,1:nj),p(1:ni,1:nj),rho(1:ni,1:nj)
double precision :: evisc(1:ni,1:nj),pnorm(1:ni,1:nj)
double precision :: uu(1:ni,1:nj),vv(1:ni,1:nj),ww(1:ni,1:nj)
double precision :: ff(1:ni,1:nj),rr(1:ni,1:nj),uv(1:ni,1:nj)
double precision :: u1l(1:ni), u2l(1:ni),uc(1:ni),uprod(1:ni,1:nj)
double precision :: momthick(1:ni),momthick0,zloc(1:ni)
double precision :: eta(1:ni,1:nj), unorm(1:ni,1:nj)
double precision :: uunorm(1:ni,1:nj),vvnorm(1:ni,1:nj),wwnorm(1:ni,1:nj)
double precision :: uvnorm(1:ni,1:nj),ffnorm(1:ni,1:nj),uumax(1:ni),temp
integer :: jlow(1:ni), jhigh(1:ni),jc(1:ni),jdelta(1:ni)
double precision :: u1n(1:ni), u2n(1:ni),pn(1:ni)
open(unit=15, file = 'statistics.dat')

read(15,*)
read(15,*)



!Now for the mixing layer domain
read(15,*)
do j=1, nj
do i=1, ni
read(15,*) x(i), y(j), &
     u(i,j),v(i,j),w(i,j),xi(i,j),p(i,j),rho(i,j), &
     evisc(i,j), uu(i,j),vv(i,j),ww(i,j),uv(i,j),ff(i,j)
enddo
enddo


!Normalisation in Browand and Latigo is done by inflow momentum thickness
!This can be computed from the velocity profile using equation 

!"Freestream" velocities are the velocities at the edges of the mixing layer
! Need a better way to do this!

do i=1,ni
do j=1, nj
if(xi(i,j) .le. 0.997) jhi = j
if(xi(i,j) .le. 0.003) jlo = j
enddo
u1l(i) = u(i,jhi)
u2l(i) = u(i,jlo)
pn(i) = p(i,jlo)
!u1l(i) = u(i,1) + 0.99*(u(i,nj) - u(i,1))
!u2l(i) = u(i,1) + 0.01*(u(i,nj) - u(i,1))
!write(*,*) i, u1l(i), u2l(i)
enddo

!open(unit=15,file='dp.dat')
!write(15,*) 'VARIABLES = "x" "u1" "u2" "dp"'
!write(15,*) 'ZONE'
!do i=1, ni
!write(15,*) x(i), u1l(i), u2l(i), (u1l(i) - u2l(i)) / (u1l(i) + u2l(i))
!enddo

!do i=1,ni
!do j=1,nj
!if(u(i,j) .le. u2l(i)) jlow(i) = j
!enddo
!do j=1,nj
!if(u(i,j) .le. u1l(i)) jhigh(i) = j
!enddo
!write(*,*) i, jlow(i), jhigh(i),u(i,jlow(i)), u(i,jhigh(i))
!enddo
!
do i=1, ni
u2n(i) = u(i,2)
u1n(i) = u(i,nj-1)
uc(i) = 0.5*(u1l(i) + u2l(i))
!write(*,*) i, u1l(i), u2l(i), uc(i)
enddo
!Compute the momentum thickness of the flow
!do i=1,ni
!do j=jlow(i),jhigh(i)
!uprod(i,j) = (u(i,jhigh(i)) - u(i,j))*(u(i,j) - u(i,jlow(i)))
!uprod(i,j) = uprod(i,j) / (u(i,jhigh(i)) - u(i,jlow(i)))**2
!write(90003,*) i,j,uprod(i,j)
!enddo
!enddo


!do i=1,ni
!do j=jlow(i),jhigh(i)
!momthick(i) = momthick(i) + 0.5*(y(j+1)-y(j))*(uprod(i,j) + uprod(i,j+1))
!enddo
!enddo

!do i=1,ni
!write(90000,*) x(i)/momthick(1), momthick(i)/momthick(2)
!enddo

!Work out locus of Uc
!do i=1,ni
!do j=1,nj
!if(u(i,j) .le. uc(i)) then
!jc(i) = j
!endif
!enddo
!zloc(i) = y(jc(i))/ momthick(1)
!enddo



!Compute self-similar flow stats
!do i=1,ni
!do j=1,nj
!eta(i,j) = (y(j) - y(jc(i))) / momthick(i)
!enddo
!enddo

do i=1,ni
do j=1,nj
unorm(i,j) = ( u(i,j)-uc(i)) / (u1l(i) - u2l(i))
uunorm(i,j) = sqrt(uu(i,j)) / (u1l(i) - u2l(i))
vvnorm(i,j) = sqrt(vv(i,j)) / (u1l(i) - u2l(i))
wwnorm(i,j) = sqrt(ww(i,j)) / (u1l(i) - u2l(i))
uvnorm(i,j) = uv(i,j) / (u1l(i) - u2l(i))**2
if(ff(i,j) .lt. 0.0) ff(i,j) = 0.0
ffnorm(i,j) = sqrt(ff(i,j))
pnorm(i,j) = (p(i,j) - pn(i)) / (1.2474*(u1l(i)-u2l(i))**2)
!write(*,*) i,j, uu(i,j), uunorm(i,j), u1l(i), u2l(i)
enddo
enddo


!do i =1, ni
!temp = 0.0
!do j=1, nj
!if(uunorm(i,j) .gt. temp) then
!temp = uunorm(i,j)
!uumax(i) = temp
!endif
!enddo
!write(90001,*) x(i), uumax(i)
!enddo 
open(unit=25, file='selfsimstats.dat')
write(25,*) 'VARIABLES = "x" "y" "u" "xi" "urms" "vrms" "wrms" "uv" "xirms" "pnorm"'
write(25,*) 'ZONE I=', ni, ', J=', nj
do j=1,nj
do i=1,ni
write(25,*) x(i),y(j), unorm(i,j), xi(i,j), uunorm(i,j),vvnorm(i,j), wwnorm(i,j), uvnorm(i,j), ffnorm(i,j),pnorm(i,j)
enddo
enddo



end program selfstats
