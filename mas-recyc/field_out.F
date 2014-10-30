!-----------------------------------------------------
! routine to output single plane and spanwise averaged
! flow visualisations
!-----------------------------------------------------
  
       subroutine field_out(x,y,z,dip,diu,djp,djv,dkp,
     &             u,v,w,sc1,rho,p,gam,step) 
  
       implicit none
! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
#include "data_type.inc"

      integer step  
      real  u(*), v(*), w(*), sc1(*), rho(*), p(*)
      real  dip(*), diu(*), djp(*),
     &      djv(*), gam(*)
      real  dkp(*), x(*),y(*), z(*)
      real  work1(maxbl2*nip1*njp1*nkp1)
      real  w2a(maxbl2*nip1*njp1), w2b(maxbl2*nip1*njp1),
     &      w2c(maxbl2*nip1*njp1), w2d(maxbl2*nip1*njp1),
     &      w2e(maxbl2*nip1*njp1), w2f(maxbl2*nip1*njp1),
     &      w2g(maxbl2*nip1*njp1), w2h(maxbl2*nip1*nkp1),
     &      w2i(maxbl2*nip1*njp1), w2j(maxbl2*nip1*nkp1)
! local variables

      real   xtemp(maxbl2*(nip1+1)*(njp1+1))
      real   ytemp(maxbl2*(nip1+1)*(njp1+1))
      integer iblk,joff,i,j,ij,itemp2(3,maxbl2)
      integer pmg_k, pmg_ij,ps_wij,itemp(3), iblk0
      integer ps_ijk,ps_wijk,p_xy,mpierr,sreq(maxbl2*3)
      integer ilev1,status(MPI_STATUS_SIZE) 
      integer ii, poff,mg_or_sg,i0,gm,fm,dm,bm,tm

      do iblk = 1,numbl

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)
         
        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1
        p_xy = mg_off_2d_xy(iblk,1) + 1
C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1 
  
       ps_wijk = 1 + ((iblk-1)*nip1*njp1*nkp1) 
       ps_wij = 1 + (iblk-1)*nip1*njp1

      call comp_schlieren(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &      rho(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &      djv(pmg_ij),dkp(pmg_k),work1(ps_wijk),1)
      call dim_ave(work1(ps_wijk),w2e,3)
	  
      call comp_schlieren(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &      rho(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &      djv(pmg_ij),dkp(pmg_k),work1(ps_wijk),2)
      call dim_ave(work1(ps_wijk),w2i,3)

       call comp_flowviz(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &      rho(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &      djv(pmg_ij),dkp(pmg_k),work1(ps_wijk),1)
       call dim_ave(work1(ps_wijk),w2f,3)   
   
   
       call comp_flowviz(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &      rho(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &      djv(pmg_ij),dkp(pmg_k),work1(ps_wijk),3)
      call dim_ave(work1(ps_wijk),w2g,3)   
   
      call dim_ave(u(ps_ijk),w2a,3)
      call dim_ave(v(ps_ijk),w2b,3)
      call dim_ave(w(ps_ijk),w2c,3)
      call dim_ave(sc1(ps_ijk),w2d,3)

      

      enddo !iblk

!      write(*,*) pe, step, 'here'

        do iblk = 1,numbl

        dm = 1000 + iblkl2g(iblk)
        fm = 3000 + iblkl2g(iblk)
        gm = 4000 + iblkl2g(iblk)
        tm = 8000 + iblkl2g(iblk)

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)

        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1
        p_xy = mg_off_2d_xy(iblk,1) + 1
C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1

       ps_wijk = 1 + ((iblk-1)*nip1*njp1*nkp1)
       ps_wij = 1 + (iblk-1)*nip1*njp1 
            call write_2dfile(w2a,w2b,w2c,w2d,
     &          w2e,w2f,w2g,w2i,x(p_xy),y(p_xy),nip1,njp1,fm,gm)

	 enddo ! iblk

 
      do iblk = 1,numbl


        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)
         
        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1
        p_xy = mg_off_2d_xy(iblk,1) + 1
C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1 
  
       ps_wijk = 1 + ((iblk-1)*nip1*njp1*nkp1) 
       ps_wij = 1 + (iblk-1)*nip1*njp1


        call comp_schlieren(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &      rho(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &      djv(pmg_ij),dkp(pmg_k),work1(ps_wijk),1)


      call dim_ave(u(ps_ijk),w2a,0)
      call dim_ave(v(ps_ijk),w2b,0)
      call dim_ave(w(ps_ijk),w2c,0)
      call dim_ave(sc1(ps_ijk),w2d,0)
      call dim_ave(work1(ps_wijk),w2e,0)
      call dim_ave(gam(ps_wijk),w2f,0)


            call write_2dfile2(w2a,w2b,w2c,w2d,
     &          w2e,w2f,x(p_xy),y(p_xy),nip1,njp1,dm)


      call dimy_ave(work1(ps_wijk), w2h,2)
      call write_xzfile(w2h,x(p_xy),z(p_xy),
     &          nip1,njp1,nkp1,tm)
      enddo


    
! write out y-z planes here	  
!--------------------------------------------------------------
      do iblk = 1,numbl

        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
        nkp1 = bl_nijkp1(3,iblk)
         
        pmg_k  = mg_off_1d_k(iblk,1)  + 1
        pmg_ij = mg_off_2d_ij(iblk,1) + 1
        p_xy = mg_off_2d_xy(iblk,1) + 1
!C POINTERS FOR SINGLE GRID QUANTITIES

        ps_ijk = sg_off(iblk) + 1 
  
       ps_wijk = 1 + ((iblk-1)*nip1*njp1*nkp1) 
       ps_wij = 1 + (iblk-1)*nip1*njp1	  
	  
        bm = 2000 + iblkl2g(iblk)
        	  
       call yz_output(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &     rho(ps_ijk),sc1(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &     djv(pmg_ij),dkp(pmg_k),x(p_xy),y(p_xy),z(p_xy),
     &     nip1,njp1,nkp1,step,bm)  

      enddo  


!write out x-z y-averaged schlieren and shadowgraph

   
       end subroutine

 
!-------------------------------------------------------------
       subroutine dim_ave(f, f2d,dim)

! amend this to fit in with the current implementation
! generalise this to a) so spanwise averaging
! b) do cross-stream averaging, and c) pick out a plane 
! at mid-span of domin

      implicit none

! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
      real f(nip1*njp1*nkp1), f2d(nip1*njp1)

      integer   ibl,mpierr
      integer   i, j, k, ijk, ij
      integer   dim
      real      sum

      if(dim .eq. 3) then
       jo = nip1
       ko = nip1*njp1
         
         do j = 1, njp1
            do i = 1, nip1
               ij = i + jo*(j-1) 
               sum = 0.0
               do k = 2, nkp1-1
                  ijk = i + jo*(j-1) + ko*(k-1) 
                  sum = sum + f(ijk)
               enddo
               f2d(ij) = sum / float(nkp1-2)
           enddo
         enddo

       else if(dim .eq. 0) then
! write out single plane from nk/2
         do j = 1, njp1
            do i = 1, nip1
               ij = i + jo*(j-1) 
                  ijk = i + jo*(j-1) + ko*((nkp1-1)/2) 
               f2d(ij) = f(ijk)        
        if(pe.eq.1) then
!        write(92,*) i,j, f2d(ij)
        endif
           enddo
         enddo
 
       endif !dim
      
	  return 
      end

!----------------------------------------------------------------
! This routine calculates schileren

      subroutine comp_schlieren(u,v,w,rho,dip,diu,djp,djv,dkp,
     &             work1,dim) 
 
      implicit none
C PARAMETERS
#include "mpif.h"
#include "dim.inc"
#include "bc_type.inc"

C INCLUDE FILES
#include "data_type.inc"
#include "cb.gridp.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "const.inc"

      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     rho(nip1*njp1*nkp1), work1(nip1*njp1*nkp1)
      real     dip(nip1,njp1), diu(nip1,njp1), djp(nip1,njp1),
     &         djv(nip1,njp1)
      real     dkp(nkp1)
!local
      integer  i, j, k, jk, ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr,
     &         ijkse, ijknw, ijkne, ijkel, ijknl, ijkwr, ijker, ijksr,
     &         ijknr,dim,ijkee,ijkrr,ijknn
      integer  koff,iside_next,bc,iblk_next,iblk,ierr
      real     drdx,drdy,drdz,red,temp
      real     d2rdx,d2rdy,d2rdz,rhoef,rhowf,rhonf,rhosf,rholf,rhorf


 
         jo = nip1
         ko = nip1*njp1
         ni = nip1 - 1
         nj = njp1 - 1
         nk = nkp1 - 1
       temp = 0.0
      do k = 2,nk
        koff = ko*(k-1)
        do j = 2,nj
          jk = jo*(j-1) + koff
          do i = 2,ni

            ijk   = i + jk

            ijkw  = ijk  - 1
            ijke  = ijk  + 1
            ijks  = ijk  - jo
            ijkn  = ijk  + jo
            ijkl  = ijk  - ko
            ijkr  = ijk  + ko

            ijkse = ijks + 1
            ijknw = ijkn - 1
            ijkne = ijkn + 1

            ijkel = ijke - ko
            ijknl = ijkn - ko
            ijkwr = ijkw + ko
            ijker = ijke + ko
            ijksr = ijks + ko
            ijknr = ijkn + ko

            rhoef = 0.5*(rho(ijk) + rho(ijke))
            rhowf = 0.5*(rho(ijkw) + rho(ijk))
            rhonf = 0.5*(rho(ijkn) + rho(ijk))
            rhosf = 0.5*(rho(ijks) + rho(ijk))
            rholf = 0.5*(rho(ijkl) + rho(ijk))
            rhorf = 0.5*(rho(ijk) + rho(ijkr))   
            drdx = (rhoef - rhowf) /diu(i,j)
            drdy = (rhonf - rhosf)/ djv(i,j)
            drdz = (rhorf - rholf)/dkp(k)

            if(dim .eq. 1) then
            work1(ijk) = sqrt(drdx**2+drdy**2+drdz**2)
            if(work1(ijk) .ge. temp) temp = work1(ijk)                     
!            write(*,*) pe, work1(ijk), temp
            endif

              enddo
            enddo
       enddo
      call MPI_REDUCE(temp, red, 1, REAL_TYPE_MPI,
     &  MPI_MAX, 0, MPI_COMM_WORLD, ierr)
!      if(master) write(*,*), red
      call MPI_BCAST(red, 1, REAL_TYPE_MPI, 0, MPI_COMM_WORLD, ierr)
!      if(master) write(*,*) red
      
      do k = 2,nk
        koff = ko*(k-1)
        do j = 2,nj
          jk = jo*(j-1) + koff
          do i = 2,ni

            ijk   = i + jk

            ijkw  = ijk  - 1
            ijke  = ijk  + 1
            ijks  = ijk  - jo
            ijkn  = ijk  + jo
            ijkl  = ijk  - ko
            ijkr  = ijk  + ko

            ijkse = ijks + 1
            ijknw = ijkn - 1
            ijkne = ijkn + 1

            ijkel = ijke - ko
            ijknl = ijkn - ko
            ijkwr = ijkw + ko
            ijker = ijke + ko
            ijksr = ijks + ko
            ijknr = ijkn + ko

            rhoef = 0.5*(rho(ijk) + rho(ijke))
            rhowf = 0.5*(rho(ijkw) + rho(ijk))
            rhonf = 0.5*(rho(ijkn) + rho(ijk))
            rhosf = 0.5*(rho(ijks) + rho(ijk))
            rholf = 0.5*(rho(ijkl) + rho(ijk))
            rhorf = 0.5*(rho(ijk) + rho(ijkr))   
            drdx = (rhoef - rhowf) /diu(i,j)
            drdy = (rhonf - rhosf)/ djv(i,j)
            drdz = (rhorf - rholf)/dkp(k)

            if(dim .eq. 1) then
            work1(ijk) = 0.8*exp(-15.d0*work1(ijk)/ red)
            endif
            if(dim .eq. 2) then
            work1(ijk) = 0.0!8*exp(-150.d0*work1(ijk)/ red)
            endif            

              enddo
            enddo
       enddo




      return 
      end subroutine 

!----------------------------------------------------------
      subroutine write_2dfile(worka,workb,workc,workd,
     &          worke,workf,workg,worki,x,y,nip1,njp1,fm,gm) 

      implicit none

      integer nip1, njp1
      real worka(*),workb(*),workc(*),workd(*),worke(*) 
	  real workf(*), workg(*), worki(*)
      real x(nip1+1,njp1+1), y(nip1+1,njp1+1)
      integer joff,j,i,ij,fm,gm

!THESE INDICES WORK FOR OUTPUTTING SPANWISE AVERAGED DATA
!MOVE BACK TO MAIN ROUTINE AND USE FOR DATA TRANSFER ETC.

      write(fm,*) 'ZONE I=', nip1-1,' J=', njp1-1
      joff = nip1
      do j=2,njp1
      do i=2, nip1
       ij = i+joff*(j-1)
      write(fm,*) x(i,j),y(i,j),worka(ij),workb(ij),
     &            workc(ij),workd(ij),worke(ij),worki(ij)
      worka(ij) = 0.0
      enddo
      enddo
	  
      write(gm,*) 'ZONE I=', nip1-1,' J=', njp1-1
      joff = nip1
      do j=2,njp1
      do i=2, nip1
       ij = i+joff*(j-1)
      write(gm,*) x(i,j),y(i,j),workf(ij),workg(ij)
      enddo
      enddo
	  
	  
	  
        return

      end subroutine
!---------------------------------------------------------

!----------------------------------------------------------
      subroutine write_2dfile2(worka,workb,workc,workd,
     &          worke,workf,x,y,nip1,njp1,dm) 

      implicit none

      integer nip1, njp1
      real worka(*),workb(*),workc(*),workd(*),worke(*) 
      real x(nip1+1,njp1+1), y(nip1+1,njp1+1), workf(*)
      integer joff,j,i,ij,dm

!THESE INDICES WORK FOR OUTPUTTING SPANWISE AVERAGED DATA
!MOVE BACK TO MAIN ROUTINE AND USE FOR DATA TRANSFER ETC.

      write(dm,*) 'ZONE I=', nip1-1,' J=', njp1-1
      joff = nip1
      do j=2,njp1
      do i=2, nip1
       ij = i+joff*(j-1)
      write(dm,*) x(i,j),y(i,j),worka(ij),workb(ij),
     &            workc(ij),workd(ij),worke(ij),workf(ij)
      enddo
      enddo
	    
	  
        return

      end subroutine
!---------------------------------------------------------
!----------------------------------------------------------
      subroutine write_xzfile(worka,
     &          x,z,nip1,njp1,nkp1,tm) 
	 
       implicit none

      integer nip1, nkp1,njp1
      real worka(*)
      real x(nip1+1,njp1+1), z(nkp1+1)
      integer joff,j,i,ik,koff,k,tm

!THESE INDICES WORK FOR OUTPUTTING SPANWISE AVERAGED DATA
!MOVE BACK TO MAIN ROUTINE AND USE FOR DATA TRANSFER ETC.
       
      write(tm,*) 'ZONE I=', nip1-2,' J=', nkp1-2
      koff = nip1
      do k=2,nkp1-1
!      write(*,*) k, z(k) 
      do i=2, nip1-1
       ik = i+koff*(k-1)
      write(tm,*) x(i,1),z(k),worka(ik)
      enddo
      enddo
	    
	  
        return

      end subroutine
!---------------------------------------------------------
!-------------------------------------------------------------
       subroutine dimy_ave(f, f2d,dim)

      implicit none

! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
      real f(nip1*njp1*nkp1), f2d(nip1*nkp1)

      integer   ibl,mpierr,koff
      integer   i, j, k, ijk, ik
      integer   dim
      real      sum

      if(dim .eq. 2) then
       jo = nip1
       koff = nip1 
       ko = nip1*njp1
         
         do k = 1, nkp1
            do i = 1, nip1
               sum = 0.0
               do j = 2, njp1-1
                  ik = i + koff*(k-1)
                  ijk = i + jo*(j-1) + ko*(k-1) 
                  sum = sum + f(ijk)
               enddo
               f2d(ik) = sum !/ float(njp1-2)
           enddo
         enddo

       endif !dim
      
	  return 
      end

!-----------------------------------------------------------

      subroutine comp_flowviz(u,v,w,rho,dip,diu,djp,djv,dkp,
     &             work1,dim) 
 
      implicit none
C PARAMETERS

#include "dim.inc"
#include "bc_type.inc"

C INCLUDE FILES

#include "cb.gridp.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "const.inc"

      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     rho(nip1*njp1*nkp1)
      real     dip(nip1,njp1), diu(nip1,njp1), djp(nip1,njp1),
     &         djv(nip1,njp1)
      real     dkp(nkp1), work1(nip1*njp1*nkp1)
!local
      integer  i, j, k, jk, ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr,
     &         ijkse, ijknw, ijkne, ijkel, ijknl, ijkwr, ijker, ijksr,
     &         ijknr,dim
      integer  koff,iside_next,bc,iblk_next,iblk
      real     dvdx,dudy,dvdz,dwdy,drdxef,drdxwf,drdysf,drdynf,
     &         drdzlf,drdzrf 


 
         jo = nip1
         ko = nip1*njp1
         ni = nip1 - 1
         nj = njp1 - 1
         nk = nkp1 - 1

      do k = 2,nk
        koff = ko*(k-1)
        do j = 2,nj
          jk = jo*(j-1) + koff
          do i = 2,ni
 
            ijk   = i + jk

            ijkw  = ijk  - 1
            ijke  = ijk  + 1
            ijks  = ijk  - jo
            ijkn  = ijk  + jo
            ijkl  = ijk  - ko
            ijkr  = ijk  + ko

            ijkse = ijks + 1
            ijknw = ijkn - 1
            ijkne = ijkn + 1

            ijkel = ijke - ko
            ijknl = ijkn - ko
            ijkwr = ijkw + ko
            ijker = ijke + ko
            ijksr = ijks + ko
            ijknr = ijkn + ko

!vorticity first
            dvdx = (v(ijk) - v(ijkw))/dip(i,j)
            dudy = (u(ijk) - u(ijks))/djp(i,j)
            dwdy = (w(ijk) - w(ijks))/djp(i,j)
            dvdz = (v(ijk) - v(ijkl))/dkp(k)

       if(dim .eq.1) work1(ijk) = dwdy-dvdz !x-vort
       if(dim .eq.3) work1(ijk) = dvdx-dudy !z-vort

 
              enddo
            enddo
       enddo

      return 
      end subroutine 
	  
!--------------------------------------------------------------------------
       subroutine yz_output(u,v,w,
     &      rho,sc1,dip,diu,djp,
     &      djv,dkp,x,y,z,
     &      nip1,njp1,nkp1,step,bm)  
 

      implicit none
#include "dim.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

      integer  i,j,k,num,ilow,nip1,njp1,nkp1,unitnum,ijk,step
      integer  joff,koff,countout,bm
      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     rho(nip1*njp1*nkp1),sc1(nip1*njp1*nkp1)
      real     dip(nip1,njp1), diu(nip1,njp1), djp(nip1,njp1),
     &         djv(nip1,njp1)
      real     dkp(nkp1)
      real     x(nip1+1,njp1+1), y(nip1+1,njp1+1), z(nkp1+1)
      logical  writeout


!work out variable based file naming here.....
!can remove the do i=2,nip1 loop
     

      countout = 0
      if((step-fiel_freq) .eq. 0) then
      do num = 1,9
      unitnum = bm
       if(yz(num) .gt. -0.2) then
       do i=2, nip1
        if(x(i,1) .lt. yz(num) .and. x(nip1,1) .gt. yz(num) 
     &     .and. x(2,1) .lt. yz(num)) then
        open(unit = unitnum)
        if(countout .eq. 0) then
        write(unitnum,*) 'VARIABLES = "z" "y" "u" "v" "w" "xi" "rho"'
        endif
        countout = countout+1
       endif
      enddo   
      endif
      enddo
      endif


!.................................................

      do num = 1,9
      writeout = .false.
      ilow = 0
       if(yz(num) .gt. -0.2) then
       do i=2, nip1
        if(x(i,1) .lt. yz(num) .and. x(nip1,1) .gt. yz(num) 
     &     .and. x(2,1) .lt. yz(num)) then
        ilow = i
        writeout = .true.
       endif
      enddo


      unitnum = bm
      joff = nip1
      if(writeout) then
      write(unitnum,*) 'ZONE I=', nkp1-2,' J=', njp1-1
      do j=2,njp1
      do k=2, nkp1-1
        koff = nip1*njp1*(k-1)
        ijk = ilow + joff*(j-1) + koff 
        write(unitnum,*) z(k),y(ilow,j),u(ijk),v(ijk),w(ijk),
     &                        sc1(ijk),rho(ijk)
      enddo
      enddo
      endif !writeout

   
      endif ! yz

      enddo  !num

       return
       end
