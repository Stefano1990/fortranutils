program maketurbbl

implicit none

integer :: i,j,k,ilow,m2,jdum,n2
parameter(m2=128,n2=m2+1)
real :: yoverdeltau(1:33), uoverU(1:33)
real :: yoverdeltaurms(1:33), urmsoverU(1:33)
real :: yoverdeltavrms(1:33), vrmsoverU(1:33)
real :: yoverdeltawrms(1:33), wrmsoverU(1:33)
real :: uinput(1:m2), urmsinput(1:m2),vrmsinput(1:m2)
real :: wrmsinput(1:m2)
real :: delta, U0, yu(1:n2), y(1:n2),temp,temp2,sum


print *, 'Type in freestream velocity'
read *, U0

print *, 'Type in target boundary layer thickness (m)'
read *, delta

!read in profiles

open (unit=10, file='lamuprofile.txt')
read(10,*)
do j=1, 33
read(10,*) yoverdeltau(j), uoverU(j)
enddo

open(unit=11,file='lamurms3percent.txt')
read(11,*)
do j=1,32
read(11,*)yoverdeltaurms(j), urmsoverU(j)
vrmsoverU(j) = urmsoverU(j) *0.6 
wrmsoveru(j) = urmsoverU(j) *0.8
urmsoverU(j) = urmsoverU(j) 
!write(*,*) urmsoverU(j), vrmsoverU(j)
enddo

! Profiles read in

!  Read in y-grid to be used in simulation
open(unit=14, file = 'upperycoord.3d')
do j=1,m2
read(14,*) jdum, y(j)
!y(j) = y(j) - 0.00004
enddo

do j=1, m2-1
yu(j) = 0.5*(y(j) + y(j+1))
enddo

! u,urms, and wrms must reside on the yu grid
! v, and vrms will reside on the y grid

!Un-normalise the y- co-ordinate for each variable
!Un-normalise the velocity variables

do j=1, 33
yoverdeltau(j) = yoverdeltau(j)*delta
yoverdeltavrms(j) = yoverdeltaurms(j)*delta
yoverdeltawrms(j) = yoverdeltaurms(j)*delta
uoverU(j)=uoverU(j)*U0
vrmsoverU(j) = vrmsoverU(j)*U0
wrmsoverU(j) = wrmsoverU(j)*U0
enddo

do j=1,32
yoverdeltaurms(j) = yoverdeltaurms(j)*delta
urmsoverU(j) = urmsoverU(j)*U0
enddo

!Interpolate data onto computational mesh
!u first

do j = 1, m2
 do i = 1, 32
  if(yoverdeltau(i) .le. yu(j)) then
   ilow = i
   endif
 enddo
   uinput(j) = uoverU(ilow) + (uoverU(ilow+1) - uoverU(ilow)) * &
              ( (yu(j) - yoverdeltau(ilow)) / (yoverdeltau(ilow+1) - yoverdeltau(ilow)))
  if(ilow .eq. 32) then
   uinput(j) = uoverU(33)
  endif
enddo

!now urms
do j = 1, m2
 do i = 1, 32
  if(yoverdeltaurms(i) .le. yu(j)) then
   ilow = i
   endif
 enddo
   urmsinput(j) =urmsoverU(ilow) + (urmsoverU(ilow+1) - urmsoverU(ilow)) * &
              ( (yu(j) - yoverdeltaurms(ilow)) / &
              (yoverdeltaurms(ilow+1) - yoverdeltaurms(ilow)))
  if(ilow .ge. 29) then
   urmsinput(j) = 0.005*U0 
  endif
enddo

!now vrms
do j = 1, m2
 do i = 1, 32
  if(yoverdeltavrms(i) .le. y(j)) then
   ilow = i
   endif
 enddo
   vrmsinput(j) = vrmsoverU(ilow) + (vrmsoverU(ilow+1) - vrmsoverU(ilow)) * &
              ( (y(j) - yoverdeltavrms(ilow)) / &
               (yoverdeltavrms(ilow+1) - yoverdeltavrms(ilow)))
  if(ilow .ge. 29) then
   vrmsinput(j) = 0.005*U0
  endif
enddo

!now wrms
do j = 1, m2
 do i = 1, 32
  if(yoverdeltawrms(i) .le. yu(j)) then
   ilow = i
   endif
 enddo
   wrmsinput(j) = wrmsoverU(ilow) + (wrmsoverU(ilow+1) - wrmsoverU(ilow)) * &
              ( (yu(j) - yoverdeltawrms(ilow)) / &
              (yoverdeltawrms(ilow+1) - yoverdeltawrms(ilow)))
  if(ilow .ge. 29) then
   wrmsinput(j) = 0.005*U0
  endif
enddo

!Write out raw velocity data to a file
! jdum, u, v, urms, vrms, wrms

open(unit = 40, file = 'inflowprofiles.3d')
do j=1, m2
if(uinput(j) .gt. U0) uinput(j) =U0 
write(40,*) j+1, uinput(j), 0.0, 0.0, urmsinput(j), vrmsinput(j), wrmsinput(j)
enddo

open(unit = 41, file = 'profiles.3d')
write(41,*) 'VARIABLES = "y" "u" "urms" "vrms" "wrms"'
do j=1, m2
if(uinput(j) .gt. U0) uinput(j) =U0
write(41,*) y(j+1), uinput(j), urmsinput(j), vrmsinput(j), wrmsinput(j)
enddo

sum = 0.0
do j=1, m2-1
!write(*,*) j, (u1jet_in(j) / u1jet_in(m2))
temp = (uinput(j) / U0) * (1.0 - uinput(j) / U0)
temp2 = (uinput(j+1) / U0) * (1.0 - uinput(j+1) / U0)
sum = sum + 0.5*(yu(j+1) - yu(j)) * (temp+temp2)
enddo

write(*,*) 'theta =' ,sum

end program



