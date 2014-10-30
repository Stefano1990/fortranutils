! outputs should be rhou,rhov,rhow, rhouu,rhovv, rhoww, rhouv
      subroutine statistics(u,v,w,p,rho,um,vm,wm,pm,rhom,ruv,rvw,
     &                  rwu,ruu,rvv,rww,
     &                  rrh,gam,gamsc1,
     &                  sc1,fm,rff,
     &                  sc1var,sc1dis,sc1varm,sc1dism,
     &                  Temperature,T_m,
     &                  dip,djp,rp,dkp,dkw)

      implicit none

#include "dim.inc"
#include "bc_type.inc"

#include "cb.lgrid.inc"
#include "cb.lgridp.inc"
#include "cb.control.inc"
      
      integer max3d, max2d

!This is rather inelegant hard-coding of an array      
      parameter (max3d = 36*66*260)
      parameter (max2d = 10000)
      
      real   u(*), v(*), w(*), p(*), gam(*), gamsc1(*)
      real   um(*), vm(*), wm(*),rhom(*),rho(*)
      real   sc1(*), fm(*),pm(*)
      real   sc1var(*), sc1dis(*), sc1varm(*), sc1dism(*)
      real   ut(max3d), vt(max3d), wt(max3d)
      real   tmp1(max3d),tmp2(max3d),tmp3(max3d)
      real   tmp4(max3d),tmp7(max3d),tmp8(max3d)
      real  w2a(maxbl2*nip1*njp1), w2b(maxbl2*nip1*njp1),
     &      w2c(maxbl2*nip1*njp1), w2d(maxbl2*nip1*njp1),
     &      w2e(maxbl2*nip1*njp1), w2f(maxbl2*nip1*njp1)

      real   rff(*),ruu(*),rvv(*),rww(*),ruv(*),rvw(*),rwu(*)
      real   dip(*),djp(*),rp(*),dkp(*),dkw(*),rrh(*)
      real   Temperature(*), T_m(*)      
      
      integer   bl_nip1(maxbl), bl_njp1(maxbl), bl_nkp1(maxbl)
      integer   i, j, k, ij, ijk, ijkw, ijke, ijks, ijkn, ijkl, ijkr
      integer   ijksw, ijkse, ijkne, ijknw, ijkwr, ijker, ijksr
      integer   ijknr, ijkel, ijknl,  ps_ijk
      integer   ibl
      integer   ioffsv, ioffxy
      
      

      do ibl = 1, numbl
         bl_nip1(ibl) = bl_nijkp1(1,ibl)
         bl_njp1(ibl) = bl_nijkp1(2,ibl)
         bl_nkp1(ibl) = bl_nijkp1(3,ibl)
      enddo

!Calculate variables at cell centre

      call centre(numbl,bl_nip1,bl_njp1,bl_nkp1,sg_off,
     &            u,v,w,ut,vt,wt)

! Compute mean statistics 
 
      do ibl = 1,numbl

        nip1 = bl_nip1(ibl)
        njp1 = bl_njp1(ibl)
        nkp1 = bl_nkp1(ibl)

        ioffsv = sg_off(ibl)
	ps_ijk = sg_off(ibl) + 1 
		
        do i = 1, nip1
         do j = 1, njp1
          do k = 1, nkp1

           ijk = ioffsv + i + nip1*(j-1) + nip1*njp1*(k-1)

       tmp1(ijk) = rho(ijk)*ut(ijk)
       tmp2(ijk) = rho(ijk)*vt(ijk)
       tmp3(ijk) = rho(ijk)*wt(ijk)
          enddo
         enddo 
        enddo
        
!spanwise average flow field
        call dim_ave(tmp1(ps_ijk),w2a,3)
        call dim_ave(tmp2(ps_ijk),w2b,3)
        call dim_ave(tmp3(ps_ijk),w2c,3)
        call dim_ave(sc1(ps_ijk),w2d,3)
        call dim_ave(p(ps_ijk),w2e,3)
        call dim_ave(rho(ps_ijk),w2f,3)


! add to statistics
!change to hace no_sam to match gentles code
          do i = 1, nip1
            do j = 1, njp1
              ij  = ioffxy + i  + nip1*(j-1)
              ijk = ioffsv + i  + nip1*(j-1) + nip1*njp1*(k-1)
              um(ij) = (um(ij)*(no_sam-1) + w2a(ij))/no_sam
              vm(ij) = (vm(ij)*(no_sam-1) + w2b(ij))/no_sam
              wm(ij) = (wm(ij)*(no_sam-1) + w2c(ij))/no_sam
              fm(ij) = (fm(ij)*(no_sam-1) + w2d(ij))/no_sam
              pm(ij) = (pm(ij)*(no_sam-1) + w2e(ij))/no_sam
              rhom(ij) = (rhom(ij)*(no_sam-1) + w2f(ij))/no_sam
          enddo
        enddo


!do spanave of sc1var, and evisc here
!T_m rather lazily used to compute mean eddy viscosity
       call dim_ave(sc1var(ps_ijk),w2a,3)
       call dim_ave(gam(ps_ijk),w2b,3)
          do i = 1, nip1
            do j = 1, njp1
              ij  = ioffxy + i  + nip1*(j-1)
              sc1varm(ij) = (sc1varm(ij)*(no_sam-1) + w2a(ij))/no_sam
              T_m(ij) = (T_m(ij)*(no_sam-1) + w2b(ij))/no_sam
!              sc1dism(ij) = (st1*(sc1dism(ij) + sc1dis(ij)) -
!     &                        st2*sc1dis(ij))*stden
            enddo
          enddo

      enddo !ibl


!Second-order statistics

      do ibl = 1,numbl
        ioffsv = sg_off(ibl)

        nip1 = bl_nip1(ibl) 
        njp1 = bl_njp1(ibl)  
        nkp1 = bl_nkp1(ibl) 
	         
        do i = 1, nip1
         do j = 1, njp1
          do k = 1, nkp1
           ijk = ioffsv + i + nip1*(j-1) + nip1*njp1*(k-1)
           tmp1(ijk) = rho(ijk)*ut(ijk)*ut(ijk)
           tmp2(ijk) = rho(ijk)*vt(ijk)*vt(ijk)
           tmp3(ijk) = rho(ijk)*wt(ijk)*wt(ijk)
           tmp4(ijk) = rho(ijk)*ut(ijk)*vt(ijk)
           tmp7(ijk) = sc1(ijk)*sc1(ijk)
           tmp8(ijk) = rho(ijk)*rho(ijk)
          enddo
         enddo
        enddo
     

        call dim_ave(tmp1(ps_ijk),w2a,3)
        call dim_ave(tmp2(ps_ijk),w2b,3)
        call dim_ave(tmp3(ps_ijk),w2c,3)
        call dim_ave(tmp4(ps_ijk),w2d,3)
        call dim_ave(tmp7(ps_ijk),w2e,3)
        call dim_ave(tmp8(ps_ijk),w2f,3)
	
!add to statistics
          do i = 1, nip1
            do j = 1, njp1
              ij  = ioffxy + i  + nip1*(j-1)
              ijk = ioffsv + i  + nip1*(j-1) + nip1*njp1*(k-1)
        ruu(ij) = (ruu(ij)*(no_sam-1) + (w2a(ij)))/no_sam
        rvv(ij) = (rvv(ij)*(no_sam-1) + (w2b(ij)))/no_sam
        rww(ij) = (rww(ij)*(no_sam-1) + (w2c(ij)))/no_sam
        ruv(ij) = (ruv(ij)*(no_sam-1) + (w2d(ij)))/no_sam 
        rff(ij) = (rff(ij)*(no_sam-1) + (w2e(ij)))/no_sam
        rrh(ij) = (rrh(ij)*(no_sam-1) + (w2f(ij)))/no_sam
          enddo
        enddo

!now for the remainder of the statistics
        do i = 1, nip1
         do j = 1, njp1
          do k = 1, nkp1
           ijk = ioffsv + i + nip1*(j-1) + nip1*njp1*(k-1)
           tmp1(ijk) = rho(ijk)*ut(ijk)*wt(ijk)
           tmp2(ijk) = rho(ijk)*vt(ijk)*wt(ijk)
          enddo
         enddo
        enddo

        call dim_ave(tmp1(ps_ijk),w2a,3)
        call dim_ave(tmp2(ps_ijk),w2b,3)

!add to statistics
          do i = 1, nip1
            do j = 1, njp1
              ij  = ioffxy + i  + nip1*(j-1)
              ijk = ioffsv + i  + nip1*(j-1) + nip1*njp1*(k-1)
        rwu(ij) = (rwu(ij)*(no_sam-1) + (w2a(ij)))/no_sam
        rvw(ij) = (rvw(ij)*(no_sam-1) + (w2b(ij)))/no_sam
          enddo
        enddo
      enddo !ibl

! Subgrid-stresses
! spanwise-averaged this routine

!      do ibl = 1,numbl
!        ioffsv = sg_off(ibl)
!        ioffxy = mg_off_2d_ij(ibl,1)

!        nip1 = bl_nip1(ibl)  
!        njp1 = bl_njp1(ibl)  
!        nkp1 = bl_nkp1(ibl) 
	         
!        do k = 2, nkp1-1
!          do i = 2, nip1-1
!            do j = 2, njp1-1

!              ijk = ioffsv + i  + nip1*(j-1) + nip1*njp1*(k-1)
!              ij  = ioffxy + i  + nip1*(j-1)

!              ijkw  = ijk  - 1
!              ijke  = ijk  + 1
!              ijks  = ijk  - nip1
!              ijkn  = ijk  + nip1
!              ijkl  = ijk  - nip1*njp1
!              ijkr  = ijk  + nip1*njp1
!              ijksw = ijks - 1
!              ijkse = ijks + 1
!              ijknw = ijkn - 1
!              ijkne = ijkn + 1
!              ijkwr = ijkw + nip1*njp1
!              ijker = ijke + nip1*njp1
!              ijksr = ijks + nip1*njp1
!              ijknr = ijkn + nip1*njp1
!              ijkel = ijke - nip1*njp1
!              ijker = ijke + nip1*njp1
!              ijknl = ijkn - ko
!              ijknr = ijkn + ko

!              usf = 0.25*(u(ijk) + u(ijke) + u(ijks) + u(ijkse))
!              unf = 0.25*(u(ijk) + u(ijke) + u(ijkn) + u(ijkne))
!              ulf = 0.25*(u(ijk) + u(ijke) + u(ijkl) + u(ijkel))
!              urf = 0.25*(u(ijk) + u(ijke) + u(ijkr) + u(ijker))

!              vwf = 0.25*(v(ijk) + v(ijkn) + v(ijkw) + v(ijknw))
!              vef = 0.25*(v(ijk) + v(ijkn) + v(ijke) + v(ijkne))
!              vlf = 0.25*(v(ijk) + v(ijkn) + v(ijkl) + v(ijknl))
!              vrf = 0.25*(v(ijk) + v(ijkn) + v(ijkr) + v(ijknr))

!              wwf = 0.25*(w(ijk) + w(ijkr) + w(ijkw) + w(ijkwr))
!              wef = 0.25*(w(ijk) + w(ijkr) + w(ijke) + w(ijker))
!              wsf = 0.25*(w(ijk) + w(ijkr) + w(ijks) + w(ijksr))
!              wnf = 0.25*(w(ijk) + w(ijkr) + w(ijkn) + w(ijknr))


!              dfdip = (fef-fwf)/dip(ij)
!              dfdjp = (fnf-fsf)/djp(ij)
!              dfdkp = (frf-flf)/(rp(ij)*dkw(k))

!              dudi = (u(ijke) - u(ijk)) / dip(ij)
!              dudj = (unf - usf) / djp(ij)
!              dudk = (urf-ulf) /  (rp(ij)*dkp(k))

!              dvdi = (vef - vwf) / dip(ij)
!              dvdj = (v(ijkn) - v(ijk)) / djp(ij)
!              dvdk = (vrf - vlf) / djp(ij)

!              dwdi = (wef-wwf) / dip(ij)
!              dwdj = (wnf-wsf) / djp(ij)
!              dwdk = (w(ijkr) - w(ijk)) / (rp(ij)*dkp(k))

!              gam11 = gam(ijk)
!              gam22 = gam(ijk)
!              gam33 = gam(ijk)
!              gam12 = gam(ijk)

!need to work out dwdi,dvdk,dudk,dwdj
!              tmp1(ijk) = -gam11*(dudi+dudi)
!              tmp2(ijk) = -gam22*(dvdj+dvdj)
!              tmp3(ijk) = -gam33*(dwdk+dwdk)
!              tmp4(ijk) = -gam12*(dudj+dvdi)
!              tmp7(ijk) = -gam11*(dwdi+dudk)
!              tmp8(ijk) = -gam11*(dwdj+dvdk)

!        call dim_ave(tmp1(ps_ijk),w2a,3)
!        call dim_ave(tmp2(ps_ijk),w2b,3)
!        call dim_ave(tmp3(ps_ijk),w2c,3)
!        call dim_ave(tmp4(ps_ijk),w2d,3)
!        call dim_ave(tmp7(ps_ijk),w2e,3)
!        call dim_ave(tmp8(ps_ijk),w2f,3)

!              suu(ij) = (suu(ij)*(no_sam-1) + (w2a(ij)))/no_sam
!              svv(ij) = (svv(ij)*(no_sam-1) + (w2b(ij)))/no_sam 
!              sww(ij) = (sww(ij)*(no_sam-1) + (w2c(ij)))/no_sam
!              suv(ij) = (suv(ij)*(no_sam-1) + (w2d(ij)))/no_sam
!              swu(ij) = (swu(ij)*(no_sam-1) + (w2e(ij)))/no_sam
!              svw(ij) = (svw(ij)*(no_sam-1) + (w2f(ij)))/no_sam
!              suu(ijk)=(st1*(suu(ijk)+tmp11(ijk))-st2*tmp11(ijk))*stden
!              svv(ijk)=(st1*(svv(ijk)+tmp22(ijk))-st2*tmp22(ijk))*stden
!              sww(ijk)=(st1*(sww(ijk)+tmp33(ijk))-st2*tmp33(ijk))*stden
!              suv(ijk)=(st1*(suv(ijk)+tmp12(ijk))-st2*tmp12(ijk))*stden
!              suf(ijk)=(st1*(suf(ijk)+tmp1f(ijk))-st2*tmp1f(ijk))*stden
!              svf(ijk)=(st1*(svf(ijk)+tmp2f(ijk))-st2*tmp2f(ijk))*stden


!              suu(ijk)=(st1*(suu(ijk)+tmp11(ijk))-st2*tmp11(ijk))*stden
!              svv(ijk)=(st1*(svv(ijk)+tmp22(ijk))-st2*tmp22(ijk))*stden
!              sww(ijk)=(st1*(sww(ijk)+tmp33(ijk))-st2*tmp33(ijk))*stden
!              suv(ijk)=(st1*(suv(ijk)+tmp12(ijk))-st2*tmp12(ijk))*stden
!              suf(ijk)=(st1*(suf(ijk)+tmp1f(ijk))-st2*tmp1f(ijk))*stden
!              svf(ijk)=(st1*(svf(ijk)+tmp2f(ijk))-st2*tmp2f(ijk))*stden

!            enddo
!          enddo
!        enddo
!      enddo !ibl


      no_sam = no_sam+1

      return

      end

      subroutine centre(numbl,bl_nip1,bl_njp1,bl_nkp1,sg_off,
     &                  u,v,w,ut,vt,wt)

      implicit none

      integer   numbl
      integer   bl_nip1(numbl), bl_njp1(numbl), bl_nkp1(numbl)
      integer   sg_off(numbl)

      real      u(*), v(*), w(*)
      real      ut(*), vt(*), wt(*)

C  LOCAL

      integer   ibl
      integer   nip1, njp1, nkp1, ni, nj, nk
      integer   jo, ko
      integer   i, j, k, ijk, ijke, ijkn, ijkr, ijks, ijkl, ijkw
      integer   ijksr, ijkse, ijkle, ijkne
      integer   ijknw, ijknr, ijknl, ijkre, ijkrw

      do ibl = 1, numbl

         nip1 = bl_nip1(ibl)
         njp1 = bl_njp1(ibl)
         nkp1 = bl_nkp1(ibl)

         ni = nip1 - 1
         nj = njp1 - 1
         nk = nkp1 - 1

         jo = nip1
         ko = nip1*njp1

         do k = 2, nk
            do j = 2, nj
               do i = 2, ni

                  ijk = i + jo*(j-1) + ko*(k-1) + sg_off(ibl)
                  ijke = ijk + 1
                  ijkn = ijk + jo
                  ijkr = ijk + ko

                  ut(ijk) = 0.5*(u(ijk) + u(ijke))
                  vt(ijk) = 0.5*(v(ijk) + v(ijkn))
                  wt(ijk) = 0.5*(w(ijk) + w(ijkr))

               enddo
            enddo
         enddo

         do k = 2, nk
            do j = 2, nj
               ijk = 1 + jo*(j-1) + ko*(k-1) + sg_off(ibl)
               ijke = ijk + 1
               ijkn = ijk + jo
               ijkne = ijkn + 1
               ijkr = ijk + ko
               ijkre = ijkr + 1
               ut(ijk) = u(ijke)
               vt(ijk) = 0.25*(v(ijk) + v(ijke) + v(ijkn) + v(ijkne))
               wt(ijk) = 0.25*(w(ijk) + w(ijke) + w(ijkr) + w(ijkre))
               ijk = nip1 + jo*(j-1) + ko*(k-1) + sg_off(ibl)
               ijkw = ijk - 1
               ijkn = ijk + jo
               ijknw = ijkn - 1
               ijkr = ijk + ko
               ijkrw = ijkr - 1
               ut(ijk) = u(ijk)
               vt(ijk) = 0.25*(v(ijk) + v(ijkw) + v(ijkn) + v(ijknw))
               wt(ijk) = 0.25*(w(ijk) + w(ijkw) + w(ijkr) + w(ijkrw))
            enddo
         enddo

         do i = 2, ni
            do j = 2, nj
               ijk = i + jo*(j-1) + ko*(1-1) + sg_off(ibl)
               ijkr = ijk + ko
               ijkn = ijk + jo
               ijknr = ijkn + ko
               ijke = ijk + 1
               ijkre = ijkr + 1
               wt(ijk) = w(ijkr)
               vt(ijk) = 0.25*(v(ijk) + v(ijkr) + v(ijkn) + v(ijknr))
               ut(ijk) = 0.25*(u(ijk) + u(ijkr) + u(ijke) + u(ijkre))
               ijk = i + jo*(j-1) + ko*(nkp1-1) + sg_off(ibl)
               ijkl = ijk - ko
               ijkn = ijk + jo
               ijknl = ijkn - ko
               ijke = ijk + 1
               ijkle = ijkl + 1
               wt(ijk) = w(ijk)
               vt(ijk) = 0.25*(v(ijk) + v(ijkl) + v(ijkn) + v(ijknl))
               ut(ijk) = 0.25*(u(ijk) + u(ijkl) + u(ijke) + u(ijkle))
            enddo
         enddo

         do k = 2, nk
            do i = 2, ni
               ijk = i + jo*(1-1) + ko*(k-1) + sg_off(ibl)
               ijkn = ijk + jo
               ijkr = ijk + ko
               ijknr = ijkr + jo
               ijke = ijk + 1
               ijkne = ijkn + 1
               vt(ijk) = v(ijkn)
               wt(ijk) = 0.25*(w(ijk) + w(ijkn) + w(ijkr) + w(ijknr))
               ut(ijk) = 0.25*(u(ijk) + u(ijkn) + u(ijke) + u(ijkne))
               ijk = i + jo*(njp1-1) + ko*(k-1) + sg_off(ibl)
               ijks = ijk - jo
               ijkr = ijk + ko
               ijksr = ijkr - jo
               ijke = ijk + 1
               ijkse = ijks + 1
               vt(ijk) = v(ijk)
               wt(ijk) = 0.25*(w(ijk) + w(ijks) + w(ijkr) + w(ijksr))
               ut(ijk) = 0.25*(u(ijk) + u(ijks) + u(ijke) + u(ijkse))
            enddo
         enddo

      enddo

      return
      end
!----------------------------------------------------------------------	  
!-------------------------------------------------
! dump stats routine added by WAM 22/06/10
! dumps spanwise-avered statistics to a file
!-------------------------------------------------

       subroutine dump_stats(x,y,z,um,vm,wm,ruv,rvw,rwu,
     &     ruu,rvv,rww,
     &     rrh,fm,rff,
     &     Temperature,T_m,pm,rhom,sc1varm)
	  
       implicit none
	   
! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
#include "data_type.inc"

       real  x(*),y(*), z(*), um(*), vm(*),wm(*),ruv(*)
       real  rvw(*),rwu(*),ruu(*),rvv(*),rww(*),rrh(*)
       real  fm(*),rff(*),Temperature(*), T_m(*),pm(*),rhom(*)
	   real  sc1varm(*)
	   
      real  w2a(maxbl2*nip1*njp1), w2b(maxbl2*nip1*njp1),
     &      w2c(maxbl2*nip1*njp1), w2d(maxbl2*nip1*njp1),
     &      w2e(maxbl2*nip1*njp1), w2f(maxbl2*nip1*njp1)
      real   xtemp(maxbl2*(nip1+1)*(njp1+1))
      real   ytemp(maxbl2*(nip1+1)*(njp1+1))
      integer iblk,joff,i,j,ij,itemp2(3,maxbl2)
      integer pmg_k, pmg_ij,ps_wij,itemp(3), iblk0
      integer ps_ijk,ps_wijk,p_xy,mpierr,sreq(maxbl2*3)
      integer ilev1,status(MPI_STATUS_SIZE) 
      integer ii, poff,mg_or_sg,i0, mm, sm
!arrays and stuff here



        do iblk = 1,numbl

         mm = 5000 + iblkl2g(iblk)
         sm = 6000 + iblkl2g(iblk)


        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
            p_xy = mg_off_2d_xy(iblk,1) + 1
 
            call write_meanfile(um,vm,wm,fm,pm,rhom,T_m,
     &          x(p_xy),y(p_xy),nip1,njp1,1,no_sam,mm)


            call write_stressfile(ruu,rvv,rww,ruv,rwu,rvw,rff,rrh,
     &      sc1varm,x(p_xy),y(p_xy),nip1,njp1,1,no_sam,sm)



        enddo !iblk

       if (master) write(*,*) 'mean and stress files written'

	  
        return
       end subroutine


!------------------------------------------------
! subroutine write_meanfile
! writes mean data to file
!------------------------------------------------


       subroutine write_meanfile(worka,workb,workc,workd,
     &          worke,workf,workg,x,y,nip1,njp1,over,no_sam,mm)
	   

      implicit none

!#include "cb.control.inc"

      integer nip1, njp1
      integer over, no_sam,mm
      real worka(*),workb(*),workc(*),workd(*),worke(*) 
	  real workf(*),workg(*)
      real x(nip1+1,njp1+1), y(nip1+1,njp1+1)
      integer joff,j,i,ij

!THESE INDICES WORK FOR OUTPUTTING SPANWISE AVERAGED DATA
!MOVE BACK TO MAIN ROUTINE AND USE FOR DATA TRANSFER ETC.
!	  if(over .eq. 1) then
          open(mm)
	  rewind(mm)
	  write(mm,*) no_sam
      write(mm,*)'VARIABLES = "x" "y" "ru" "rv" "rw" "sc1" "p" "r" "Ev"'
!	  endif

      write(mm,*) 'ZONE I=', nip1-2,', J=', njp1-2
      joff = nip1
      do j=2,njp1-1
      do i=2, nip1-1
       ij = i+joff*(j-1)
      write(mm,*) x(i,j),y(i,j),worka(ij),workb(ij),
     &            workc(ij),workd(ij),worke(ij),workf(ij),workg(ij)
      enddo
      enddo
  
	  
       close(mm)
        return

      end subroutine
	   
!------------------------------------------------
! subroutine write_stressfile
! writes mean data to file
!------------------------------------------------


       subroutine write_stressfile(worka,workb,workc,workd,
     &          worke,workf,workg,workh,worki,x,y,nip1,njp1,
     &          over,no_sam,sm)
	   
	
      implicit none

      integer nip1, njp1
      integer over, no_sam,sm
      real worka(*),workb(*),workc(*),workd(*),worke(*) 
	  real workf(*),workg(*),workh(*),worki(*)
      real x(nip1+1,njp1+1), y(nip1+1,njp1+1)
      integer joff,j,i,ij



         open(sm)
	  rewind(sm)
	  write(sm,*) no_sam
      write(sm,*)'VARIABLES = "x" "y" "ruu" "rvv" "rww" "ruv" 
     & "rwu" "rvw" "ff" "rr" "sc1var"'


      write(sm,*) 'ZONE I=', nip1-2,', J=', njp1-2
      joff = nip1
      do j=2,njp1-1
      do i=2, nip1-1
       ij = i+joff*(j-1)
      write(sm,*) x(i,j),y(i,j),worka(ij),workb(ij),
     &            workc(ij),workd(ij),worke(ij),workf(ij),workg(ij),
     &            workh(ij),worki(ij)
      enddo
      enddo
  
       close(sm)
        return

      end subroutine
  
!-------------------------------------------------------------
! subroutine spectra
! writes out instantaneous flow variables along 
! mixing layer centreline
!-------------------------------------------------------------

        subroutine spectra(x,y,z,dip,diu,djp,djv,dkp,
     &             u,v,w,sc1,rho)
	 
       implicit none
! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
#include "data_type.inc"
  
      real  u(*), v(*), w(*), sc1(*), rho(*)
      real  dip(*), diu(*), djp(*),
     &      djv(*)
      real  dkp(*), x(*),y(*), z(*)
	  
       integer iblk,joff,i,j,ij
      integer pmg_k, pmg_ij,ps_wij,itemp(3), iblk0
      integer ps_ijk,ps_wijk,p_xy
      integer ii, poff,mg_or_sg,unitnum
	  
	  
      do iblk = 1,numbl
!generalise this?  You know that it'll always be
!the plane of the splitter plate
        if(iblkl2g(iblk) .ge. 88 .and. iblkl2g(iblk) .le. 112) then

        unitnum = 7000+iblkl2g(iblk)
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
	  
        ps_ijk = sg_off(iblk) + 1 
  
       ps_wijk = 1 + ((iblk-1)*nip1*njp1*nkp1) 
       ps_wij = 1 + (iblk-1)*nip1*njp1	  
	  
       call writefluc(u(ps_ijk),v(ps_ijk),w(ps_ijk),
     &     rho(ps_ijk),sc1(ps_ijk),dip(pmg_ij),diu(pmg_ij),djp(pmg_ij),
     &     djv(pmg_ij),dkp(pmg_k),x(p_xy),y(p_xy),z(p_xy),
     &     nip1,njp1,nkp1,unitnum)  

       endif
      enddo  
	  
		  
        return
       end subroutine
	   
!-------------------------------------------------------
! subroutien writefluc
! writes out fluctuations along cenreline in each block
!-------------------------------------------------------

      subroutine writefluc(u,v,w,
     &      rho,sc1,dip,diu,djp,
     &      djv,dkp,x,y,z,
     &      nip1,njp1,nkp1,unitnum)  
 

      implicit none
#include "dim.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

      integer  i,j,k,num,ilow,nip1,njp1,nkp1,unitnum,ijk,jk
      integer  joff,koff
      real     u(nip1*njp1*nkp1), v(nip1*njp1*nkp1), w(nip1*njp1*nkp1)
      real     rho(nip1*njp1*nkp1),sc1(nip1*njp1*nkp1)
      real     dip(nip1,njp1), diu(nip1,njp1), djp(nip1,njp1),
     &         djv(nip1,njp1)
      real     dkp(nkp1)
      real     x(nip1+1,njp1+1), y(nip1+1,njp1+1), z(nkp1+1)
      logical  writeout
  
 
 

        koff = nip1*njp1*((nkp1/2)-1) 
        jk = nip1*(3-1) + koff
         do i=2, nip1-1
         ijk = i + jk
         write(unitnum,*) i, u(ijk), v(ijk), w(ijk), sc1(ijk), rho(ijk)
         enddo

  
        return
  
       end subroutine
       
!--------------------------------------------------------------
!subroutine read stats
!read the 2-d stats from old files to continue accumulation
!wam 13/09/13
!--------------------------------------------------------------
!----------------------------------------------------------------------   


       subroutine read_stats(um,vm,wm,ruv,rvw,rwu,
     &     ruu,rvv,rww,
     &     rrh,fm,rff,
     &     Temperature,T_m,pm,rhom,sc1varm)
        
       implicit none
        
! include files here
#include "dim.inc"
#include "cb.lgridp.inc"
#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "comm.inc"
#include "cb.comm.inc"
#include "mpif.h"
#include "data_type.inc"

       real  um(*), vm(*),wm(*),ruv(*)
       real  rvw(*),rwu(*),ruu(*),rvv(*),rww(*),rrh(*)
       real  fm(*),rff(*),Temperature(*), T_m(*),pm(*),rhom(*)
       real sc1varm(*)

      integer iblk,joff,i,j,ij,itemp2(3,maxbl2)
      integer pmg_k, pmg_ij,ps_wij,itemp(3), iblk0
      integer ps_ijk,ps_wijk,p_xy,mpierr,sreq(maxbl2*3)
      integer ilev1,status(MPI_STATUS_SIZE)
      integer ii, poff,mg_or_sg,i0, mm, sm
!arrays and stuff here



        do iblk = 1,numbl

         mm = 5000 + iblkl2g(iblk)
         sm = 6000 + iblkl2g(iblk)

!         write(*,*) pe, mm
        nip1 = bl_nijkp1(1,iblk)
        njp1 = bl_nijkp1(2,iblk)
            p_xy = mg_off_2d_xy(iblk,1) + 1       
       

            call read_meanfile(um,vm,wm,fm,pm,rhom,T_m,
     &          nip1,njp1,1,no_sam,mm)


            call read_stressfile(ruu,rvv,rww,ruv,rwu,rvw,rff,rrh,
     &      sc1varm,nip1,njp1,1,no_sam,sm)



        enddo !iblk

       if (master) write(*,*) 'mean and stress files read in'

        
        return
       end subroutine


!------------------------------------------------
! subroutine read_meanfile
! reads mean data from file
!------------------------------------------------


       subroutine read_meanfile(worka,workb,workc,workd,
     &          worke,workf,workg,nip1,njp1,over,no_sam,mm)
        

      implicit none

!#include "cb.control.inc"

      integer nip1, njp1
      integer over, no_sam,mm
      real worka(*),workb(*),workc(*),workd(*),worke(*)
      real workf(*), workg(*)
      real dum(nip1+1,njp1+1)
      integer joff,j,i,ij

!THESE INDICES WORK FOR OUTPUTTING SPANWISE AVERAGED DATA
!MOVE BACK TO MAIN ROUTINE AND USE FOR DATA TRANSFER ETC.
          open(mm)
          rewind(mm)
          read(mm,*) no_sam
      read(mm,*)


      read(mm,*)
      joff = nip1
      do j=2,njp1-1
      do i=2, nip1-1
       ij = i+joff*(j-1)
      read(mm,*) dum(i,j),dum(i,j),worka(ij),workb(ij),
     &            workc(ij),workd(ij),worke(ij),workf(ij),workg(ij)
      enddo
      enddo

        
       close(mm)
        return

      end subroutine
        
!------------------------------------------------
! subroutine read_stressfile
! reads stress data from file
!------------------------------------------------

       subroutine read_stressfile(worka,workb,workc,workd,
     &          worke,workf,workg,workh,worki,nip1,njp1,over,no_sam,sm)
        
        
      implicit none

      integer nip1, njp1
      integer over, no_sam,sm
      real worka(*),workb(*),workc(*),workd(*),worke(*)
          real workf(*),workg(*),workh(*),worki(*)
      real dum(nip1+1,njp1+1)
      integer joff,j,i,ij

!THESE INDICES WORK FOR OUTPUTTING SPANWISE AVERAGED DATA
!MOVE BACK TO MAIN ROUTINE AND USE FOR DATA TRANSFER ETC.
         open(sm)
          rewind(sm)
          read(sm,*)
      read(sm,*)


      read(sm,*)
      joff = nip1
      do j=2,njp1-1
      do i=2, nip1-1
       ij = i+joff*(j-1)
      read(sm,*) dum(i,j),dum(i,j),worka(ij),workb(ij),
     &            workc(ij),workd(ij),worke(ij),workf(ij),workg(ij),
     &             workh(ij),worki(ij)
      enddo
      enddo

       close(sm)
        return

      end subroutine       
       

