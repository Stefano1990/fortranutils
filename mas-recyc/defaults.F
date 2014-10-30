       SUBROUTINE Defaults
!      -------------------
         implicit none

#include "dim.inc"
#include "cb.control.inc"
#include "cb.comm.inc"
#include "cb.gridp.inc"
!include cb.file here
           
         integer i, ierr
         
         
      ierr = 0

      if (rotate.ne.0) then
        gotrad = .true.
      else
        gotrad = .false.
      end if

 
         step_no = 1
         visc = 1.8E-05
         rhohi = 1.24
         rholo = 1.24
         csles = 0.1
         last = 10
         dtime = 0.000001
         New_Run = .TRUE.
         numvar = 47 
         t_relax = 0.00001
         ncycle(1) = 1
         ncycle(2) = 1
         ncycle(3) = 1
         ncycle(4) = 20
         ncycle(5) = 1
         numsolve = 5 
         Stat_Freq = 10**4
         Fiel_Freq = 10**4
         Dump_Freq = 10**4      
         Spec_Freq = 10**4
         no_sam = 1
         statstart = 500000
         addback = 1
         refpresblock = 1
         jsplit = 1000
         iturb = 1
         ca    = 2.3
         ck    = 0.41
         aplus = 25.0
         inlet1 = 1
         inlet2 = 2
         rhoini = 1.2474
         debugging = .false.
         specout = .false.
         writewall = .false.
         psolver=101
         relax_pmg = 1.6
         relax_pmgc = 1.6
         tolc = 5.0e-11
         tolf = 5.0e-13
         numpre = 1
         numpost = 1
	 numcits = 20
	 nummgits = 5
	 nummgl = 0
	 meth_csolve = 0
         meth_smooth = 1
	 meth_cycle = 1
         ststep = 1
         relax(1) = 1
         relax(2) = 1
         relax(3) = 1
         relax(4) = 1
         relax(5) = 1  
	 do i=1, maxvar
	 ifile(i) = 0
	 enddo

         do i=1,9
         yz(i) = -1.0
         enddo

         ifskip = 100000
         
      inbcglength   = 0
      numinbcg      = 0
      inbcgblk      = 1
      inbcgplane    = 3
      inbcgfileskip = 0
      inbcgfileoff  = 0
      do i = 1,maxvar
        inbcgvar(i) = 0
      end do
         RETURN
 
       END
