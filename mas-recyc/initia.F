C=======================================================================
C Set initial condition
C=======================================================================

       subroutine initia(u,v,w,sc1,p,rho,rhonm1,rhonp1,y,uprof
     &                    ,uuprof,vvprof,wwprof)

      implicit none

#include "dim.inc"
#include "const.inc"
#include "les.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "comm.inc"
#include "cb.comm.inc"

C ARGUMENTS
      real         u(*),p(*),rho(*),rhonm1(*),rhonp1(*)
      real         v(*)
      real         w(*), uprof(*),uuprof(*),vvprof(*),wwprof(*)
      real         sc1(*)
      real         y(*)

C LOCAL
      integer iseed
      real rhoin
      integer nip1, njp1, nkp1, ni, nj, nk, ibl
      integer i, j, k, ij, ijk, jo, ko
      integer boff1, boff3
      integer joff, koff, joffy, jkoff,dm

      real    perturb, random, gasdev
      real    uin1, vin1, win1, uin2, vin2, win2, uin, vin, win, utotal
      real    sigu1, sigu2, sigu, fin1, fin2, fin
      real    R0, Acoef

C-----------------------------------------------------------------------
C INITIALIZATION
C-----------------------------------------------------------------------

      fin1 = 1.0
      fin2 = 0.0


      do ibl = 1, numbl

         nip1 = mg_nijkp1(1,ibl,1)
         njp1 = mg_nijkp1(2,ibl,1)
         nkp1 = mg_nijkp1(3,ibl,1)

         ni = nip1 - 1
         nj = njp1 - 1
         nk = nkp1 - 1

         jo = nip1
         ko = nip1*njp1

         boff1 = mg_off_2d_xy(ibl,1)
         boff3 = sg_off(ibl)

         do k = 1, nkp1
            koff = (k-1)*ko + boff3

            do j = 1, njp1
               joff = (j-1)*jo
               jkoff= joff + koff
            joffy = joff + boff1

               do i = 1, nip1
                  ijk = i + jkoff
              ij = i + joffy

       if(iblkl2g(ibl) .le. 24) then
       call random_number(random)
       perturb = 0.05*uprof(2)*2.0*(random-0.5)
        u(ijk) = uprof(j) + perturb
       perturb = 0.05*uprof(2)*2.0*(random-0.5)
        v(ijk) = perturb
       perturb = 0.05*uprof(2)*2.0*(random-0.5)
        w(ijk) = perturb
        sc1(ijk) = 0.0
        rho(ijk) = rholo
        rhonm1(ijk) = rholo
        rhonp1(ijk) = rholo
       else if(iblkl2g(ibl) .ge. 25 .and. iblkl2g(ibl) .le. 48) then
        call random_number(random)
       perturb = 0.05*uprof(2)*2.0*(random-0.5)
        u(ijk) = uprof(j+(nj-1)) + perturb
       perturb = 0.05*uprof(2)*2.0*(random-0.5)
        v(ijk) = perturb
       perturb = 0.05*uprof(2)*2.0*(random-0.5)
        w(ijk) = perturb
        sc1(ijk) = 0.0
        rho(ijk) = rholo
        rhonm1(ijk) = rholo
        rhonp1(ijk) = rholo 
        else if(iblkl2g(ibl) .ge. 49 .and. iblkl2g(ibl) .le. 72) then 
        call random_number(random)
       perturb = 0.1*uprof(2*(nj-2))*2.0*(random-0.5)
!         write(*,*) pe, perturb
        u(ijk) = uprof(j) + perturb
       perturb = 0.1*uprof(2*(nj-2))*2.0*(random-0.5)
        v(ijk) = perturb
       perturb = 0.1*uprof(2*(nj-2))*2.0*(random-0.5)
        w(ijk) = perturb
        sc1(ijk) = 1.0
        rho(ijk) = rhohi
        rhonm1(ijk) = rhohi
        rhonp1(ijk) = rhohi  
        else if(iblkl2g(ibl) .ge. 73 .and. iblkl2g(ibl) .le. 96) then 
        call random_number(random)
       perturb = 0.1*uprof(2*(nj-2))*2.0*(random-0.5)
        u(ijk) = uprof(j+(nj-1)) + perturb
       perturb = 0.1*uprof(2*(nj-2))*2.0*(random-0.5)
        v(ijk) = perturb
       perturb = 0.1*uprof(2*(nj-2))*2.0*(random-0.5)
        w(ijk) = perturb
        sc1(ijk) = 1.0
        rho(ijk) = rhohi
        rhonm1(ijk) = rhohi
        rhonp1(ijk) = rhohi
       endif                   
               enddo
            enddo
         enddo

      end do

      if(debugging) then
      write(*,*) pe, '-Flow field initialised-'
      endif
      return
      end


      function random(s)
c
c     function to create "random" sequence of numbers between 0 and 0.1
c
c     input:   s      - number between 0 and 1
c
c     output:  random - number between 0 and 1
c              s      - s=rando
c
      save
      random=100.0e0*exp(s)
      random=random-float(int(random))
      s=random
      return
      end

      FUNCTION GASDEV(IDUM)
C
C  THIS FUNTION RETRUNS A NORMAL DISTRIBUTION WITH ZERO MEAN AND
C  UNIT VARIANCE
C  USES RAN1(IDUM) AS A SOURCE OF UNIFORM DEVIATES
C  IDUM HAS TO BE A NEGATIVE INTEGER TO INITIALISE
C
      INTEGER IDUM
      REAL GASDEV
      INTEGER ISET
      REAL FAC,GSET,RSQ,V1,V2,RAN1
      SAVE ISET,GSET
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
    1    V1=2.*RAN1(IDUM)-1.
         V2=2.*RAN1(IDUM)-1.
         RSQ=V1**2+V2**2
         IF (RSQ.GE.1..OR.RSQ.EQ.0.) GOTO 1
         FAC=SQRT(-2.*LOG(RSQ)/RSQ)
         GSET=V1*FAC
         GASDEV=V2*FAC
         ISET=1
      ELSE
         GASDEV=GSET
         ISET=0
      END IF
C
      RETURN
      END
C
      FUNCTION RAN1(IDUM)
      INTEGER IDUM,IA,IM,IQ,IR,NTAB,NDIV
      REAL RAN1,AM,EPS,RNMX
      PARAMETER (IA=16807, IM=2147483647, AM=1./IM, IQ=127773, IR=2836,
     1           NTAB=32,NDIV=1+(IM-1)-NTAB, EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER J,K,IV(NTAB),IY
      SAVE IV,IY
      DATA IV /NTAB*0/, IY /0/
C
      IF (IDUM.LE.0 .OR. IY.EQ.0) THEN
          IDUM=MAX(-IDUM,1)
          DO J=NTAB+8,1,-1
            K=IDUM/IQ
            IDUM=IA*(IDUM-K*IQ)-IR*K
            IF (IDUM.LT.0) IDUM=IDUM+IM
            IF (J.LE.NTAB) IV(J)=IDUM
          ENDDO
          IY=IV(1)
      END IF
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF (IDUM.LT.0) IDUM=IDUM+IM
      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=MIN(AM*IY,RNMX)
C
      RETURN
      END
C
