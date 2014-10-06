C=======================================================================
C LOOP THROUGH BLOCKS TO CALCULATE SOURCE TERM FOR PRESSURE EQUATION
C=======================================================================


      subroutine sour_p(diu,djv,rdiv,rdju,didjp,dkp,rdidjp,dkw, 
     &  u,v,w,rho,rhonp1,Drhonp1Dt,s,srhot,srhoux,gia,gja)

      implicit none

#include "mpif.h"
#include "dim.inc"
#include "data_type.inc"
#include "comm.inc"

#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

      real         diu(*), djv(*), rdiv(*), rdju(*), didjp(*)
      real         dkp(*)

      real         rdidjp(*) 
      real         dkw(*)   

      real         u(*), v(*), w(*), rho(*),rhonp1(*)  
      real         s(*), srhot(*), srhoux(*)

      real         Drhonp1Dt(*)  
      real         gia(*), gja(*)

C LOCAL

      integer      iblk
      integer      nip1, njp1, nkp1
      integer      io, jo, ko
      integer      p1, p2, p3, p3s
      real         rdum(2)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      do 10 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

        io = 1
        jo = nip1
        ko = nip1*njp1

C POINTERS

        p1  = mg_off_1d_k(iblk,1)  + 1
        p2  = mg_off_2d_ij(iblk,1) + 1
        p3  = mg_off(iblk,1)       + 1
        p3s = sg_off(iblk)         + 1


        call sour_p_block(nip1,njp1,nkp1,io,jo,ko,dtime,
     &       diu(p2),djv(p2),rdiv(p2),rdju(p2),didjp(p2),dkp(p1),
C     &       u(p3s),v(p3s),w(p3s),rho(p3),s(p3),
     &       u(p3s),v(p3s),w(p3s),rho(p3),rhonp1(p3),Drhonp1Dt(p3),
     &       s(p3),srhot(p3),srhoux(p3),gia,gja,
     &       rdidjp(p2), dkw(p1),iblk) 
 10   continue

C!!!! debug

c     call file_3(20,'s',1,s,1,1,rdum(1),rdum(2))

      return
      end
