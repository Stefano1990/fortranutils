PROGRAM PDF

INCLUDE 'yzjet3d.inc'

INTEGER :: i,j,k,iloc, t, tmax, dpsi, etaj1, etaj2, dim, stream
REAL ::  phi7max, scalar(1:M1,1:M2,0:1000), ximy(1:M1,1:M2), Cm(1:M1,1:M2)
REAL :: Cmy(1:M1,1:M2), y05(1:M1), eta(1:M1,1:M2)
real ::  rl(1:m1), theta1, prod(1:m1,1:m2,0:1000), total(1:m1,1:m2), pm(1:m1,1:m2)
real :: ximm(1:m1),pmy(1:m1),ev(1:m1),px(1:m1),x0,normy(1:m1,1:m2)
integer :: code, norm, hestream
integer :: jlow(1:m1), jhigh(1:m1), jdiff(1:m1)
real :: yupper(1:m1), ylower(1:m1), molehe, moleair
real :: uc(1:m1), rho1, rho2, Rel(1:m1),rhoparam,pmyy(1:m1)

!Read in data from lean reactant in high-speed stream simulation






PRINT *, 'Enter number of zones in input file.'

READ *, tmax


WRITE(*,*) 'Prcoessing....'

dpsi=25


open(unit=50, file = 'yzflowviz.dat')
read(50,*)


phi7max = 1000.0
!read in flow output at each time instant and put values in scalar bin

do t=1, tmax

read(50,*)
do j=1, m2
do i=1,m1
read(50,*) x(i), y(j), u1(i,j,1),u2(i,j,1), u3(i,j,1), pass1(i,j,1)
!pass1_acc(i,j) = pass1_acc(i,j) + pass1(i,j,1)
pass1(i,j,1) = pass1(i,j,1) * 1000.0

enddo
enddo



DO I=1, M1
DO J=1, M2
DO K = 0, 1000, dpsi
IF (pass1(i,j,1) .GE. real(K) .AND. pass1(i,j,1) .LT. real(K+dpsi)) THEN
scalar(i,j,k) = scalar(i,j,k) + 1.0
ENDIF
ENDDO
ENDDO
ENDDO


if(mod(t, 100) .eq.0)  then
WRITE(*,*) 'Completed output', t
endif

enddo

WRITE(*,*) 'averaging'

do j = 1, m2-1
yu(j) = 0.5*(y(j) + y(j+1))
enddo
yu(m2) = y(m2-1) + (yu(m2-1) - yu(m2-2))

do i=1, m1
do j=1, m2
do k=0, 1000,dpsi
scalar(i,j,k) = scalar(i,j,k) / real(tmax)
enddo
enddo
enddo

do i=1, m1
do j=1, m2
do k=0, 1000, dpsi
total(i,j) = total(i,j) + real(dpsi)/1000.0 * scalar(i,j,k)
enddo
enddo
enddo


do i=1,m1
do j=1, m2
do k=0, 1000, dpsi
scalar(i,j,k) = scalar(i,j,k) / total(i,j)
enddo
enddo
enddo





OPEN(UNIT=35, FILE='pdfplot.dat')
WRITE(35,*) 'VARIABLES = "y (m)" "<greek>f" "p"'
DO i=1, m1
WRITE(35,*) 'ZONE I=', (1000/dpsi)+1, ', J=',m2 
DO J=1,m2
DO K = 0, 1000,dpsi
if(scalar(i,j,k) .ge. 8.0) scalar(i,j,k) = 8.0
WRITE(35,*) y(j), real(K)/1000.0, scalar(i,j,k)
ENDDO
ENDDO
enddo



STOP

END PROGRAM
