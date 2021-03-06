!--------------------------------
! Impose inlet boundary condition
!--------------------------------

      subroutine fixed_inlet(hu,hv,hw,u_int,v_int,w_int,
     &        u,v,w,p,sc1,rho,gam,iseed,
     &        y,uprof,vprof,wprof,uuprof,vvprof,wwprof,umean,
     &        vmean,wmean,uumean,vvmean,wwmean,step)



      implicit none

#include 'mpif.h'
#include 'data_type.inc'

#include "dim.inc"
#include "les.inc"
#include "comm.inc"

#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc"

c Arguments

      integer iseed,step

      real    hu(*), hv(*), hw(*)
      real    u_int(*), v_int(*), w_int(*)
      real    u(*), v(*), w(*), sc1(*)
      real    p(*), rho(*), gam(*),wprof(*)
      real    y(*), uprof(*), vprof(*),uuprof(*),vvprof(*),wwprof(*)
      real    umean(*),vmean(*),wmean(*),uumean(*),vvmean(*),wwmean(*)
      
c Local

      integer iva, ibl
      integer nip1, njp1, nkp1
      integer ni, nj, nk
      integer j, k, ij, ijk, ijkw,i,joffy
      integer jo, ko,koff,jk,ijk2,jkoff
      integer boff1, boff3
      real    perturb, random, gasdev,cpfilt
      real    uin1, vin1, win1, uin2, vin2, win2, uin, vin, win, utotal
      real    sigu1, sigu2, sigu, fin1, fin2, fin
      real    rhoin, meanu, meanv, meanw,meanuu,meanvv,meanww
      real    R0, Acoef,small
      integer i1, joff, dm,jdum
      integer status(MPI_STATUS_SIZE), request
      integer ierr


! Read inlet data
      small = 0.0000001
      fin1 = 1.0
      fin2 = 0.0
      cpfilt = 0.01
! should make these blocks parameters in the input file
! for generalisation purposes
      do ibl = 1, numbl
        if ((iblkl2g(ibl) .eq. inlet1)
     &    .or.(iblkl2g(ibl) .eq. inlet2)
     &    .or.(iblkl2g(ibl).eq. inlet3)
     &    .or.(iblkl2g(ibl) .eq. inlet4)) then

            nip1 = mg_nijkp1(1,ibl,1)
            njp1 = mg_nijkp1(2,ibl,1)
            nkp1 = mg_nijkp1(3,ibl,1)

            ni = nip1 - 1
            nj = njp1 - 1
            nk = nkp1 - 1

            jo = nip1
            ko = nip1*njp1

            boff1 = mg_off_2d_xy(ibl,1)
	   i1= 2
            boff3 = sg_off(ibl)

         do j = 1, njp1
            joff = (j-1)*nip1
            ij= i1+ joff + boff1

            do k = 1, nkp1
            if (iblkl2g(ibl) == inlet1) then
                uin = uprof(j)
                vin = vprof(j)
                call random_number(random)
                perturb = 0.0025*uprof(3)*2.0*(random-0.5)
                fin = 1.0
                rhoin = rholo
                call random_number(random)
                perturb = 0.0025*uprof(3)*2.0*(random-0.5)
                win = perturb
            endif
            if (iblkl2g(ibl) == inlet2) then
                uin = uprof(j+(nj-1))
                vin = vprof(j+(nj-1))
                call random_number(random)
                perturb = 0.0025*uprof(3)*2.0*(random-0.5)
                fin = 1.0
                rhoin = rholo
                call random_number(random)
                perturb = 0.0025*uprof(3)*2.0*(random-0.5)
                win = perturb
            endif
            if (iblkl2g(ibl) == inlet3) then
                uin = uprof(j+2*(nj-1))
                vin = vprof(j+2*(nj-1))
                call random_number(random)
                perturb = 0.0025*uprof(4*(nj-2))*2.0*(random-0.5)
                fin = 0.0
                rhoin = rhohi
                call random_number(random)
                perturb = 0.0025*uprof(4*(nj-2))*2.0*(random-0.5)
                win = perturb
            endif
            if (iblkl2g(ibl) == inlet4) then
                uin = uprof(j+3*(nj-1))
                vin = vprof(j+3*(nj-1))
                call random_number(random)
                perturb = 0.0025*uprof(4*(nj-2))*2.0*(random-0.5)
                fin = 0.0
                rhoin = rhohi
                call random_number(random)
                perturb = 0.0025*uprof(4*(nj-2))*2.0*(random-0.5)
                win = perturb
            endif
             
!              if(j .lt. jsplit) then
!                uin = uprof(j)
!                vin = vprof(j)
!                call random_number(random)
!                perturb = 0.005*uprof(3)*2.0*(random-0.5)
!                fin = 0.0
!                rhoin = rholo
!                call random_number(random)
!                perturb = 0.005*uprof(3)*2.0*(random-0.5)
!                win = perturb
!              else 
!                uin = uprof(j)
!                vin = vprof(j)
!                call random_number(random)
!                perturb = 0.005*uprof(nj-1)*2.0*(random-0.5) 
!                fin = 1.0 
!                rhoin = rhohi
!                call random_number(random)
!                perturb = 0.005*uprof(nj-1)*2.0*(random-0.5)
!                win = perturb
!              endif
                  uin = uin + perturb
                  vin = vin + perturb
           ! endif
               ijk = 2 + jo*(j-1) + ko*(k-1) + boff3
               ijkw = ijk - 1
               u(ijk) = uin
               v(ijkw) = vin 
               w(ijkw) = win 

               sc1(ijkw) = fin
	       rho(ijkw) = rhoin

               u_int(ijk) = u(ijk)
               v_int(ijkw) = v(ijkw)
               w_int(ijkw) = w(ijkw)
               enddo
            enddo


        end if


      enddo
       
        
      return
      end

!---------------------------------------------------
! subroutine profile read
! reads a file containing velocity profile 
! information for inlet condition
!---------------------------------------------------

       subroutine read_profile(uprof, vprof,wprof,
     &            uuprof,vvprof,wwprof)
  
       implicit none
#include "dim.inc"
#include "les.inc"
#include "comm.inc"

#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc" 
!array defs
  
      real uprof(*), vprof(*),uuprof(*),vvprof(*),wwprof(*)
      real wprof(*)
!local
      integer ibl, j,njp1,nj,jdum,m
    
  
 !this can be changed to reflect the blocks in the inflow file above.
 !needs to be generalised  

      do ibl = 1, numbl
  
  
        m = iblkl2g(ibl)
        njp1 = mg_nijkp1(2,ibl,1)
        nj = njp1 - 1
        open(unit=10000+m,file='inflowbc.3d')
        do j = 1,(4*nj)-2
          read(10000+m,*) jdum, uprof(j), vprof(j)
        enddo
        uprof(1) = uprof(2)
        vprof(1) = vprof(2)
        close(unit=10000+m) 

      enddo
 
      return
   
      end subroutine
