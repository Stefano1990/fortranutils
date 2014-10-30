C=======================================================================
C Impose inlet boundary conditions
C=======================================================================

      subroutine recycling(hu,hv,hw,u_int,v_int,w_int,
     &        u,v,w,p,sc1,rho,gam,iseed,
     &        y,uprof,vprof,wprof,uuprof,vvprof,wwprof,umean,
     &        vmean,wmean,uumean,vvmean,wwmean,step,inblock)



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

      integer iseed,step,inblock

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
      real    utemp(1:66*258), vtemp(1:66*258), wtemp(1:66*258)
      real    uwork(1:66,1:258), vwork(1:66,1:258), wwork(1:66,1:258)
      real    perturb, random, gasdev,cpfilt
      real    uin1, vin1, win1, uin2, vin2, win2, uin, vin, win, utotal
      real    sigu1, sigu2, sigu, fin1, fin2, fin
      real    rhoin, meanu, meanv, meanw,meanuu,meanvv,meanww
      real    R0, Acoef,small
      integer i1, joff, dm,jdum
      integer status(MPI_STATUS_SIZE), request
      integer ierr,jmult


c Read inlet data
      small = 0.0000001
      fin1 = 1.0
      fin2 = 0.0
      cpfilt = 0.144
      
!      if(inblock .eq. inlet3) jmult = 0
!      if(inblock .eq. inlet4) jmult = njp1-2      
!      if(inblock .eq. inlet4) write(*,*) njp1-2, nj
! should make these blocks parameters in the input file
! for generalisation purposes
 
       if(step .eq. 1 .and. New_Run .eq. .false.) then
         do ibl = 1, numbl
         if(iblkl2g(ibl).ge.inblock.and.iblkl2g(ibl).le.inblock+23) then
         dm=iblkl2g(ibl)
         open(unit=15000+dm)
          rewind(15000+dm)
           do j=1, njp1
            read(15000+dm,*) jdum, umean(j), vmean(j), wmean(j)
            read(15000+dm,*) jdum, uumean(j), vvmean(j), wwmean(j)
           enddo
          close(15000+dm)
         endif  !iblk
        enddo !ibl
        endif !step 
             


         do ibl = 1, numbl
         if(iblkl2g(ibl).ge.inblock.and.iblkl2g(ibl).le.inblock+23) then
            dm = iblkl2g(ibl)
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
!Do avereaging part of recycling and rescaling method - step 1 in Wang and McGuirk
!      write(*,*) 'doing averaging at sampling plane'
        jk = 1
       do k = 2, nkp1
       koff = (k-1)*ko + boff3
        do j = 2, njp1  
        joff = (j-1)*jo
        jkoff= joff + koff
        joffy = joff + boff1
        meanu = 0.0
        meanv = 0.0
        meanw = 0.0
         do i = 2, nip1!  upper limit is the sampling plane
         ijk = i + jkoff
         meanu = meanu + u(ijk)
         meanv = meanv + v(ijk)
         meanw = meanw + w(ijk)
         enddo ! i   
         uwork(j,k) = meanu / real(nip1-1)
         vwork(j,k) = meanv / real(nip1-1)
         wwork(j,k) = meanw / real(nip1-1)
!mean velocities
        enddo
       enddo
  
  
        do j= 2, njp1
        meanu = 0.0
        meanv = 0.0
        meanw = 0.0   
         do k=2, nkp1
          meanu = meanu + uwork(j,k)
          meanv = meanv + vwork(j,k)
          meanw = meanw + wwork(j,k)
        enddo
        meanu = meanu / float(nkp1-1)
        meanv = meanv / float(nkp1-1)
        meanw = meanw / float(nkp1-1)                
        umean(j) = cpfilt*meanu + (1. - cpfilt)*umean(j)
        vmean(j) = cpfilt*meanv + (1. - cpfilt)*vmean(j)
        wmean(j) = cpfilt*meanw + (1. - cpfilt)*wmean(j)       
        enddo
                
         
!rms velocities
        jk = 1
       do k = 2, nkp1
       koff = (k-1)*ko + boff3
        do j = 2, njp1  
        joff = (j-1)*jo
        jkoff= joff + koff
        joffy = joff + boff1
        meanuu = 0.0
        meanvv = 0.0
        meanww = 0.0
         do i=2, nip1
         ijk = i+jkoff
         meanuu = meanuu + (u(ijk) - umean(j))**2
         meanvv = meanvv + (v(ijk) - vmean(j))**2
         meanww = meanww + (w(ijk) - wmean(j))**2         
         enddo

         uwork(j,k) = meanuu / real(nip1-1)
         vwork(j,k) = meanvv / real(nip1-1)
         wwork(j,k) = meanww / real(nip1-1)
        enddo
       enddo
       
        do j= 2, njp1
        meanuu = 0.0
        meanvv = 0.0
        meanww = 0.0   
         do k=2, nkp1
          meanuu = meanuu + uwork(j,k)
          meanvv = meanvv + vwork(j,k)
          meanww = meanww + wwork(j,k)
        enddo
        meanuu = meanuu / float(nkp1-1)
        meanvv = meanvv / float(nkp1-1)
        meanww = meanww / float(nkp1-1)
         uumean(j) = sqrt(cpfilt*meanuu 
     &         + (1.-cpfilt)*uumean(j)**2)
         vvmean(j) = sqrt(cpfilt*meanvv
     &         + (1.-cpfilt)*vvmean(j)**2)
         wwmean(j) = sqrt(cpfilt*meanww 
     &         + (1.-cpfilt)*wwmean(j)**2)
       
        enddo       

!       write(*,*) 'done with averaging'

        if(mod(step, 1000) .eq. 0) then
        if(iblkl2g(ibl).ge.inblock.and.iblkl2g(ibl).le.inblock+23) then
        dm = iblkl2g(ibl)
        open(unit=15000+dm)
        rewind(15000+dm)
        do j=1, njp1
        write(15000+dm,*) j, umean(j), vmean(j), wmean(j)
        write(15000+dm,*) j, uumean(j), vvmean(j), wwmean(j)
        enddo
        close(15000+dm)
        endif !iblkl2g
        endif  !mod
        
       endif !iblkl2g
       enddo ! iblk

!rescale the block up to the sampling plane

      
!        call MPI_BCAST(umean, 66, REAL_TYPE_MPI, inblock+7-1, 
!     &  MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(vmean, 66, REAL_TYPE_MPI, inblock+7-1, 
!     &  MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(wmean, 66, REAL_TYPE_MPI, inblock+7-1, 
!     &  MPI_COMM_WORLD, ierr)      
!         call MPI_BCAST(uumean, 66, REAL_TYPE_MPI, inblock+7-1, 
!     &  MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(vvmean, 66, REAL_TYPE_MPI, inblock+7-1, 
!     &  MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(wwmean, 66, REAL_TYPE_MPI, inblock+7-1, 
!     &  MPI_COMM_WORLD, ierr)      
      
!This loop needs to be modified so that the j planes correspond to the
!correct part of the inflow file in the upper block
!wam 13/1/14
       
! Rescale all virtual blocks       
      do ibl = 1, numbl
       if(mod(step, 500) .eq. 0 .and. step .ge. 200) then
          if (iblkl2g(ibl) .ge. inblock .and. iblkl2g(ibl) .le.
     &              inblock+23)  then
            dm = iblkl2g(ibl)
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
      if(inblock .eq. inlet3) jmult = 0
      if(inblock .eq. inlet2) jmult = njp1-2
!      if(inblock .eq. inlet4) write(*,*) njp1-2, nj
       do k = 2, nkp1
       koff = (k-1)*ko + boff3
        do j = 2, njp1
        joff = (j-1)*jo
        jkoff= joff + koff
        joffy = joff + boff1
         if(iblkl2g(ibl) .ne. inblock) then
         do i = 1, nip1!  change njp1 to end where the recycling plane is
         ijk = i + jkoff
        u(ijk) = ( uuprof(j+jmult)/(uumean(j)+small)) *
     &          (u(ijk) - umean(j)) + uprof(j+jmult)
        v(ijk) = ( vvprof(j+jmult)/(vvmean(j)+small)) *
     &          (v(ijk) - vmean(j)) + vprof(j+jmult)
        w(ijk) = ( wwprof(j+jmult)/(wwmean(j)+small)) *
     &          (w(ijk) - wmean(j)) + wprof(j+jmult)
         enddo
         else if(iblkl2g(ibl) .eq. inblock) then
          do i = 2, nip1!  change njp1 to end where the recycling plane is
         ijk = i + jkoff
        u(ijk) = ( uuprof(j+jmult)/(uumean(j)+small)) *
     &          (u(ijk) - umean(j)) + uprof(j+jmult)
        v(ijk) = ( vvprof(j+jmult)/(vvmean(j)+small)) *
     &          (v(ijk) - vmean(j)) + vprof(j+jmult)
        w(ijk) = ( wwprof(j+jmult)/(wwmean(j)+small)) *
     &          (w(ijk) - wmean(j)) + wprof(j+jmult)
         enddo  
         endif !iblkl2g for loops      
        enddo
        enddo
!        write(*,*) 'done with rescaling'
        endif  !iblkl2g
       endif !mod

!this is for the recycling of data onto the inlet plane
          if (iblkl2g(ibl) .eq. inblock+23)  then
           dm = iblkl2g(ibl)
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
         jk=1
       do k = 2, nkp1
       koff = (k-1)*ko + boff3
        do j = 2, njp1
        joff = (j-1)*jo
        jkoff= joff + koff
        joffy = joff + boff1
!jk counter is for arrays containing the mean information       
         ijk = 3 + jkoff
         ijkw = ijk-1                       
        utemp(jk) = u(ijk)
        vtemp(jk) = v(ijkw)
        wtemp(jk) = w(ijkw)
        jk = jk+1
       enddo
      enddo
      endif  !inblock+22
      enddo   !ibl - master loop

!Broadcast data to all cpus (lazy as inlet block is the target)
        call MPI_BCAST(utemp, 66*258, REAL_TYPE_MPI, inblock+23-1, 
     &  MPI_COMM_WORLD, ierr)
        call MPI_BCAST(vtemp, 66*258, REAL_TYPE_MPI, inblock+23-1, 
     &  MPI_COMM_WORLD, ierr)
        call MPI_BCAST(wtemp, 66*258, REAL_TYPE_MPI, inblock+23-1, 
     &  MPI_COMM_WORLD, ierr)


! Recycle velocity at sampling plane onto inlet plane

         do ibl = 1, numbl
          if (iblkl2g(ibl) .eq. inblock)  then         
           dm = iblkl2g(ibl)
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
         jk=1
       do k = 2, nkp1
       koff = (k-1)*ko + boff3
        do j = 2, njp1
        joff = (j-1)*jo
        jkoff= joff + koff
        joffy = joff + boff1
!jk counter is for arrays containing the mean information
         ijk = 2 + jkoff                
               ijkw = ijk - 1
                u(ijk) = utemp(jk)
                v(ijkw) = vtemp(jk)
                w(ijkw) = wtemp(jk)
               sc1(ijkw) = 1.0
              rho(ijkw) = rhohi              
               u_int(ijk) = u(ijk)
               v_int(ijkw) = v(ijkw)
               w_int(ijkw) = w(ijkw)
                jk=jk+1
               enddo
            enddo
!        write(*,*) 'done with recycling'
          endif !iblk = 1
        enddo  !ibl

       
        
      return
      end

