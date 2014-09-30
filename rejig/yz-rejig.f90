program flowrearrange


implicit none
  integer :: bl1,t,tmax,i,j,ni,nj,bl,blnm
  parameter (ni = 152, nj = 33) ! ni = Nz, nj = Ny
  double precision :: z(1:ni,1:nj,1:4),y(1:ni,1:nj,1:4)
  double precision :: u(1:ni,1:nj,1:4),v(1:ni,1:nj,1:4),w(1:ni,1:nj,1:4)
  double precision :: xi(1:ni,1:nj,1:4),rho(1:ni,1:nj,1:4)
  double precision :: zp(1:ni,1:nj*4),yp(1:ni,1:nj*4)
  double precision :: up(1:ni,1:nj*4),vp(1:ni,1:nj*4),wp(1:ni,1:nj*4)
  double precision :: xip(1:ni,1:nj*4),rhop(1:ni,1:nj*4) 
  
open(unit=16, file='yz-flowviz.dat')
write(16,*) 'VARIABLES = "z" "y" "u" "v" "w" "xi" "rho"'

print *, 'first block?'
read *,bl1

print *, 'Number of frames?'
read *,tmax

t = 1
bl = 0
  
do i=0,3
 write(*,*) 2000+bl1+i*16
 open(unit=2000+bl1+i*16)
 read(2000+bl1+i*16,*) ! Skip variables
enddo

do t=1,tmax
print *, 'Time: ',t
do bl=1,4 ! to increase by 16 every time.
 blnm = 2000+bl1+((bl-1)*16) ! Block nr.
 read(blnm,*)

 do j=1,nj
  do i=1,ni
   read(blnm,*) z(i,j,bl),y(i,j,bl),u(i,j,bl),v(i,j,bl),w(i,j,bl),xi(i,j,bl),rho(i,j,bl)
  enddo
 enddo
enddo

! make 2D arrays.
! first row.
do j=1,nj
 do i=1,ni
  zp(i,j) = z(i,j,1)
  yp(i,j) = y(i,j,1)
  up(i,j) = u(i,j,1)
  vp(i,j) = v(i,j,1)
  wp(i,j) = w(i,j,1)
  xip(i,j) = xi(i,j,1)
  rhop(i,j) = rho(i,j,1)
 enddo
enddo
! Second row.
do j=1,nj
 do i=1,ni
  zp(i,j+nj) = z(i,j,2)
  yp(i,j+nj) = y(i,j,2)
  up(i,j+nj) = u(i,j,2)
  vp(i,j+nj) = v(i,j,2)
  wp(i,j+nj) = w(i,j,2)
  xip(i,j+nj) = xi(i,j,2)
  rhop(i,j+nj) = rho(i,j,2)
 enddo
enddo
! Third row.
do j=1,nj
 do i=1,ni
  zp(i,j+(2*nj)) = z(i,j,3)
  yp(i,j+(2*nj)) = y(i,j,3)
  up(i,j+(2*nj)) = u(i,j,3)
  vp(i,j+(2*nj)) = v(i,j,3)
  wp(i,j+(2*nj)) = w(i,j,3)
  xip(i,j+(2*nj)) = xi(i,j,3)
  rhop(i,j+(2*nj)) = rho(i,j,3)
 enddo
enddo
! Forth row.
do j=1,nj
 do i=1,ni
  zp(i,j+(3*nj)) = z(i,j,4)
  yp(i,j+(3*nj)) = y(i,j,4)
  up(i,j+(3*nj)) = u(i,j,4)
  vp(i,j+(3*nj)) = v(i,j,4)
  wp(i,j+(3*nj)) = w(i,j,4)
  xip(i,j+(3*nj)) = xi(i,j,4)
  rhop(i,j+(3*nj)) = rho(i,j,4)
 enddo
enddo

! Write 2D arrays.
write(16,*) 'ZONE I=',ni,'J=',nj*4
do j=1,nj*4
 do i=1,ni
  write(16,*) zp(i,j),yp(i,j),up(i,j),vp(i,j),wp(i,j),xip(i,j),rhop(i,j)
 enddo
enddo
enddo ! time
endprogram
