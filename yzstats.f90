program yzmakestats

implicit none

integer :: i,j,k,nj,nk,t,tmax
parameter(nj = 256, nk =256)
integer :: jhi(1:nk), jlow(1:nk)
double precision :: z(1:nk), y(1:nj)
double precision :: u(1:nj,1:nk),v(1:nj,1:nk),w(1:nj,1:nk),sc1(1:nj,1:nk)
double precision :: umean(1:nj,1:nk),vmean(1:nj,1:nk),wmean(1:nj,1:nk),sc1mean(1:nj,1:nk)
double precision :: ufluc(1:nj,1:nk),vfluc(1:nj,1:nk),wfluc(1:nj,1:nk)
double precision :: uufluc(1:nj,1:nk),vvfluc(1:nj,1:nk),wwfluc(1:nj,1:nk)
double precision :: uvfluc(1:nj,1:nk),uwfluc(1:nj,1:nk),vwfluc(1:nj,1:nk)
double precision :: sc1fluc(1:nj,1:nk),uc(1:nk),u001(1:nk),u099(1:nk)
double precision :: uumean(1:nj,1:nk),vvmean(1:nj,1:nk),wwmean(1:nj,1:nk)
double precision :: uvmean(1:nj,1:nk),uwmean(1:nj,1:nk),vwmean(1:nj,1:nk)
double precision :: y001(1:nk), y099(1:nk),yc(1:nk)
double precision :: delta(1:nk),temp, U0,temp1,temp2,U01,U99


open(unit=15, file='yzflowviz.dat')
open(unit=16, file='fullyzflowviz.dat')
open(unit=17, file='timetrace.dat')
write(17,*) 'VARIABLES = "t" "y" "z" "u" "v" "w" "sc1"'
open(unit=20,file='centreline.dat')
open(unit=21,file='thickness.dat')
read(15,*)
write(16,*) 'VARIABLES = "z" "y" "u" "v" "w" "uu" "vv" "ww" "uv" "uw" "vw"'


print *, 'number of zones?'
read *, tmax
write(17,*) 'ZONE I=', nk,', J=', nj,', K=', tmax


do t= 1, tmax
! work out mean information
read(15,*) !zone header

do j=1, nj
do k=1, nk
read(15,*) z(k), y(j), u(j,k), v(j,k), w(j,k), sc1(j,k)
umean(j,k) = umean(j,k) + u(j,k)
vmean(j,k) = vmean(j,k) + v(j,k)
wmean(j,k) = wmean(j,k) + w(j,k)
sc1mean(j,k) = sc1mean(j,k) + sc1(j,k)
if(t .eq. tmax) then
umean(j,k) = umean(j,k) / real(tmax)
vmean(j,k) = vmean(j,k) / real(tmax)
wmean(j,K) = wmean(j,K) / real(tmax)
sc1mean(j,k) = sc1mean(j,k) / real(tmax)
endif
enddo
enddo

enddo !t

temp = 0.0
temp1 = 0.0
temp2 = 0.0
do k=1, nk
uc(k) = 0.5*(umean(13,k) + umean(243,k))
temp = temp + uc(k)
u001(k) = umean(13,k) + 0.01*(umean(243,k) - umean(13,k)) 
u099(k) = umean(13,k) + 0.99*(umean(243,k) - umean(13,k))
temp1 = temp1 + u001(k)
temp2 = temp2+ u099(k)
enddo

U0 = temp / real(nk) !mean centreline velocity
U01 = temp1 / real(nk) !mean low-speed vel
U99 = temp2 / real(nk) !mean hig-speed vel
write(*,*) U0, U01, U99


write(*,*) 'done with intial read in'

!work out mixing layer thickness across span, and centreline across span


!centreline first
do k=1, nk
do j=3, nj -3
if(umean(j,k) .le. U0) then
jlow(k) = j
endif
enddo !j
!write(*,*) k,jlow(k)
!linear interpolate to where value lies
yc(k) = y(jlow(k)) + (( y(jlow(k)+1) - y(jlow(k))) * (( U0 -umean(jlow(k),k)) &
           / (umean(jlow(k)+1,k) - umean(jlow(k),k))))


write(20,*) z(k), yc(k)
do j=3, nj -3
if(umean(j,k) .le. U01) then
jlow(k) = j
endif
if(umean(j,k) .le. U99) then
jhi(k) = j
endif
enddo! j
!linear interpolate to find value position
y001(k) = y(jlow(k)) + (( y(jlow(k)+1) - y(jlow(k))) * (( U01 -umean(jlow(k),k)) &
           / (umean(jlow(k)+1,k) - umean(jlow(k),k))))

y099(k) = y(jhi(k)) + (( y(jhi(k)+1) - y(jhi(k))) * (( U99 - umean(jhi(k),k)) &
           / (umean(jhi(k)+1,k) - umean(jhi(k),k))))

!y05(i) =  y(jlow) + ( ( y(jlow+1) - y(jlow)) *( (0.5 - phi7_acc(i,jlow))&
!                 / ( phi7_acc(i,jlow+1) - phi7_acc(i,jlow))))

delta(k) = ABS(y001(k)) + ABS(y099(k))
write(21,*) z(k), delta(k)
enddo !k

!write(20,*) 'ZONE'
!do k=1, nk
!write(20,*) z(k), y001(k)
!enddo
!write(20,*) 'ZONE'
!do k=1, nk
!write(20,*) z(k), y099(k)
!enddo
!now get fluctuation data

rewind(15)
read(15,*)

do t=1, tmax
read(15,*) !header info

do j=1, nj
do k=1, nk 
read(15,*) z(k), y(j), u(j,k), v(j,k), w(j,k), sc1(j,k)
ufluc(j,k) = u(j,k) - umean(j,k)
vfluc(j,k) = v(j,k) - vmean(j,k)
wfluc(j,k) = w(j,k) - wmean(j,k)
sc1fluc(j,k) = sc1(j,k) - sc1mean(j,k)
enddo
enddo

do j=1,nj
do k=1, nk
uufluc(j,k) = ufluc(j,k)**2
vvfluc(j,k) = vfluc(j,k)**2
wwfluc(j,k) = wfluc(j,k)**2
uvfluc(j,k) = ufluc(j,k)*vfluc(j,k)
uwfluc(j,k) = ufluc(j,k)*wfluc(j,k)
vwfluc(j,k) = vfluc(j,k)*wfluc(j,k)
uumean(j,k) = uumean(j,k) + uufluc(j,k)
vvmean(j,k) = vvmean(j,k) + vvfluc(j,k)
wwmean(j,k) = wwmean(j,k) + wwfluc(j,k)
uvmean(j,k) = uvmean(j,k) + uvfluc(j,k)
uwmean(j,k) = uwmean(j,k) + uwfluc(j,k)
vwmean(j,k) = vwmean(j,k) + vwfluc(j,k)
enddo
enddo

!write out instanteous fluctuation data
WRITE(16,*) 'ZONE I=', nj,', J=', nk
do j=1, nj
do k=1, nk
write(16,*) z(k), y(j), u(j,k), v(j,k), w(j,k), uufluc(j,k), &
            vvfluc(j,k), wwfluc(j,k), uvfluc(j,k), uwfluc(j,k), &
            vwfluc(j,k)
enddo
enddo


!write out pseudo-3d trace data

do j=1,nj
do k=1,nk
write(17,*) 0.0000006*1000*(real(t-1)),y(j), z(k), u(j,k), v(j,k), w(j,k), sc1(j,k)
enddo
enddo

enddo !t

!write out mean data

WRITE(16,*) 'ZONE I=', nj,', J=', nk
do j=1, nj
do k=1, nk
uumean(j,k) = uumean(j,k) / (real(tmax))
vvmean(j,k) = vvmean(j,k) / (real(tmax))
wwmean(j,k) = wwmean(j,k) / (real(tmax))
uvmean(j,k) = uvmean(j,k) / (real(tmax))
uwmean(j,k) = uwmean(j,k) / (real(tmax))
vwmean(j,k) = vwmean(j,k) / (real(tmax))
write(16,*) z(k), y(j), umean(j,k), vmean(j,k), wmean(j,k), uumean(j,k), &
            vvmean(j,k), wwmean(j,k), uvmean(j,k), uwmean(j,k), &
            vwmean(j,k)
enddo
enddo








end program
