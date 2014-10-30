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
      real    utemp(1:64*258), vtemp(1:64*258), wtemp(1:64*258)
      real    uwork(1:64,1:258), vwork(1:64,1:258), wwork(1:64,1:258)
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
      cpfilt = 0.01
      
      if(inblock .eq. inlet3) jmult = 0
      if(inblock .eq. inlet4) jmult = njp1-2      
! should make these blocks parameters in the input file
! for generalisation purposes
 
      do ibl = 1, numbl          
       if (iblkl2g(ibl) .eq. inlet3 .or. iblkl2g(ibl) .eq. inlet4) then         
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
       do k = 2, nkp1
       koff = (k-1)*ko + boff3
        do j = 2, njp1
        joff = (j-1)*jo
        jkoff= joff + koff
        joffy = joff + boff1
!jk counter is for arrays containing the mean information
         ijk = 2 + jkoff
               ijkw = ijk - 1
        call random_number(random)
       perturb = 0.05*uprof((2*nj)-2)*2.0*(random-0.5)        
        if(iblkl2g(ibl) .eq. inlet3) then 
                u(ijk) = uprof(j) + perturb
                v(ijkw) = vprof(j) + perturb
        else if(iblkl2g(ibl) .eq. inlet4) then
                u(ijk) = uprof(j+(nj-1)) + perturb
                v(ijkw) = vprof(j+(nj-1)) + perturb
        endif   
                w(ijkw) = perturb
               sc1(ijkw) = 0.0
              rho(ijkw) = rholo              
               u_int(ijk) = u(ijk)
               v_int(ijkw) = v(ijkw)
               w_int(ijkw) = w(ijkw)
               enddo
            enddo
!        write(*,*) 'done with recycling'
          endif !iblk = 1 or 33 
        enddo  !ibl


       
        
      return
      end

