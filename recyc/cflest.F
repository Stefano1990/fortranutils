C=======================================================================
C CALCULATE THE LARGEST CFL AND CFL NUMBERS
C=======================================================================

C      subroutine cflest(u,v,w,sc1,dip,djp,dkw,rp,gam)
      subroutine cflest(u,v,w,sc1,rho,Temperature,dip,djp,dkw,rp,gam)  

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "bc_type.inc"
#include "comm.inc"
#include "const.inc"
#include "data_type.inc"

C COMMON BLOCKS
#include "cb.gridp.inc"
#include "cb.lgrid.inc"
#include "cb.control.inc"
#include "cb.lcontrol.inc"
#include "cb.comm.inc"


C INPUT/OUTPUT PARAMETERS

      real    u(*), v(*), w(*), gam(*), sc1(*)
      real    rho(*), Temperature(*)  
      real    dip(*), djp(*), dkw(*), rp(*)

C LOCAL QUANTITIES

      integer iblk
      integer istr, jstr, kstr, nip1, njp1, nkp1
      integer p1, p2, p3
      integer imxcfl, jmxcfl, kmxcfl, imxdfs, jmxdfs, kmxdfs
      integer bmxcfl, bmxdfs
      integer imaxsc1, jmaxsc1, kmaxsc1
      integer iminsc1, jminsc1, kminsc1

      integer iblkmax, iblkmin

      integer imaxrho, jmaxrho, kmaxrho, iminrho, jminrho, kminrho

      integer imaxT, jmaxT, kmaxT, iminT, jminT, kminT

      integer iblkmaxrho, iblkminrho, iblkmaxT, iblkminT
      real    cfl_bl, dfs_bl, cflmax, dfsmax
      real    u_min, v_min, w_min,  u_max, v_max, w_max
      real    sc1_min, sc1_max
      real    sc1_min_bl, sc1_max_bl

      real    rho_min, rho_max, T_min, T_max

      real    rho_min_bl,rho_max_bl, T_min_bl, T_max_bl

      integer ijk_mxcfl_bl(3), ijk_mxdfs_bl(3)
      integer ijk_maxsc1_bl(3), ijk_minsc1_bl(3)

      integer ijk_maxrho_bl(3),ijk_minrho_bl(3) 

      integer ijk_maxT_bl(3),ijk_minT_bl(3)  
      real    uvw_max_bl(3), uvw_min_bl(3)

      integer ierr, ired(4), iblk0
      real    red
      integer status(MPI_STATUS_SIZE), sreq

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      cflmax = -large_number
      dfsmax = -large_number
      u_max  = -large_number
      v_max  = -large_number
      w_max  = -large_number
      sc1_max = -large_number

      rho_max = -large_number 

      T_max = - large_number 
      u_min  =  large_number
      v_min  =  large_number
      w_min  =  large_number
      sc1_min = large_number

      rho_min = large_number 

      T_min = large_number   
      iblkmax = 0
      iblkmin = 0
      iblkmaxrho = 0
      iblkminrho = 0
      iblkmaxT = 0
      iblkminT = 0

      do 5 iblk = 1,numbl

        nip1 = mg_nijkp1(1,iblk,1)
        njp1 = mg_nijkp1(2,iblk,1)
        nkp1 = mg_nijkp1(3,iblk,1)

        istr = nijkse(1,iblk,1)
        jstr = nijkse(3,iblk,1)
        kstr = nijkse(5,iblk,1)

        p1 = mg_off_1d_k(iblk,1) + 1
        p2 = mg_off_2d_ij(iblk,1) + 1
        p3 = mg_off(iblk,1) + 1

        call cflest_block(u(p3),v(p3),w(p3),sc1(p3),
     &  rho(p3),Temperature(p3), 
     &    dtime,visc,dip(p2),
     &    djp(p2),rp(p2),dkw(p1),gam(p3),istr,jstr,kstr,nip1,njp1,nkp1,
     &       cfl_bl,dfs_bl,ijk_mxcfl_bl,ijk_mxdfs_bl,uvw_min_bl,
     &       uvw_max_bl,sc1_min_bl, sc1_max_bl,
     &  rho_min_bl,rho_max_bl,T_min_bl, T_max_bl, 
     &       ijk_maxsc1_bl,ijk_minsc1_bl,
     &  ijk_maxrho_bl,ijk_minrho_bl,ijk_maxT_bl,ijk_minT_bl) 

        if (cflmax .lt. cfl_bl) then
          cflmax = cfl_bl
          imxcfl = ijk_mxcfl_bl(1)
          jmxcfl = ijk_mxcfl_bl(2)
          kmxcfl = ijk_mxcfl_bl(3)
          bmxcfl = iblkl2g(iblk)
        end if

        if (dfsmax .lt. dfs_bl) then
          dfsmax = dfs_bl
          imxdfs = ijk_mxdfs_bl(1)
          jmxdfs = ijk_mxdfs_bl(2)
          kmxdfs = ijk_mxdfs_bl(3)
          bmxdfs = iblkl2g(iblk)
        end if

        if (sc1_max .lt. sc1_max_bl) then
          sc1_max = sc1_max_bl
          imaxsc1 = ijk_maxsc1_bl(1)
          jmaxsc1 = ijk_maxsc1_bl(2)
          kmaxsc1 = ijk_maxsc1_bl(3)
          iblkmax = iblkl2g(iblk)
        endif

        if (sc1_min .gt. sc1_min_bl) then
          sc1_min = sc1_min_bl
          iminsc1 = ijk_minsc1_bl(1)
          jminsc1 = ijk_minsc1_bl(2)
          kminsc1 = ijk_minsc1_bl(3)
          iblkmin = iblkl2g(iblk)
        endif

        if (rho_max .lt. rho_max_bl) then
          rho_max = rho_max_bl
          imaxrho = ijk_maxrho_bl(1)
          jmaxrho = ijk_maxrho_bl(2)
          kmaxrho = ijk_maxrho_bl(3)
          iblkmaxrho = iblkl2g(iblk)
        endif

        if (rho_min .gt. rho_min_bl) then
          rho_min = rho_min_bl
          iminrho = ijk_minrho_bl(1)
          jminrho = ijk_minrho_bl(2)
          kminrho = ijk_minrho_bl(3)
          iblkminrho = iblkl2g(iblk)
        endif

        if (T_max .lt. T_max_bl) then
          T_max = T_max_bl
          imaxT = ijk_maxT_bl(1)
          jmaxT = ijk_maxT_bl(2)
          kmaxT = ijk_maxT_bl(3)
          iblkmaxT = iblkl2g(iblk)
        endif

        if (T_min .gt. T_min_bl) then
          T_min = T_min_bl
          iminT = ijk_minT_bl(1)
          jminT = ijk_minT_bl(2)
          kminT = ijk_minT_bl(3)
          iblkminT = iblkl2g(iblk)
        endif

        u_max  = max(u_max,uvw_max_bl(1))
        v_max  = max(v_max,uvw_max_bl(2))
        w_max  = max(w_max,uvw_max_bl(3))
        sc1_max  = max(sc1_max,sc1_max_bl)

        rho_max  = max(rho_max,rho_max_bl)

        T_max  = max(T_max,T_max_bl)
        u_min  = min(u_min,uvw_min_bl(1))
        v_min  = min(v_min,uvw_min_bl(2))
        w_min  = min(w_min,uvw_min_bl(3))
        sc1_min  = min(sc1_min,sc1_min_bl)
        rho_min  = min(rho_min,rho_min_bl)
        T_min  = min(T_min,T_min_bl)
 5    continue

      red = cflmax
      call MPI_ALLREDUCE(red, cflmax, 1, REAL_TYPE_MPI,
     &  MPI_MAX, MPI_COMM_WORLD, ierr)

      iblk0 = numbl_total + 1
      if(red .eq. cflmax) iblk0 = bmxcfl
      ired(4) = iblk0
      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,
     &  MPI_MIN, MPI_COMM_WORLD, ierr)
      if(bl2pe(iblk0).eq.pe)then
        ired(1) = imxcfl
        ired(2) = jmxcfl
        ired(3) = kmxcfl
        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,
     &    1, MPI_COMM_WORLD, sreq, ierr)
      endif
      if(pe.eq.0)then
        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &    1, MPI_COMM_WORLD, status, ierr)
        imxcfl = ired(1)
        jmxcfl = ired(2)
        kmxcfl = ired(3)
        bmxcfl = ired(4)
      endif
      if(bl2pe(iblk0).eq.pe)then
        call MPI_WAIT(sreq,status,ierr)
      endif

      red = dfsmax
      call MPI_ALLREDUCE(red, dfsmax, 1, REAL_TYPE_MPI,
     &  MPI_MAX, MPI_COMM_WORLD, ierr)

      iblk0 = numbl_total + 1
      if(red .eq. dfsmax) iblk0 = bmxdfs
      ired(4) = iblk0
      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,
     &  MPI_MIN, MPI_COMM_WORLD, ierr)
      if(bl2pe(iblk0).eq.pe)then
        ired(1) = imxdfs
        ired(2) = jmxdfs
        ired(3) = kmxdfs
        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,
     &    1, MPI_COMM_WORLD, sreq, ierr)
      endif
      if(pe.eq.0)then
        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &    1, MPI_COMM_WORLD, status, ierr)
        imxdfs = ired(1)
        jmxdfs = ired(2)
        kmxdfs = ired(3)
        bmxdfs = ired(4)
      endif
      if(bl2pe(iblk0).eq.pe)then
        call MPI_WAIT(sreq,status,ierr)
      endif

      red = sc1_max
      call MPI_ALLREDUCE(red, sc1_max, 1, REAL_TYPE_MPI,
     &  MPI_MAX, MPI_COMM_WORLD, ierr)

      iblk0 = numbl_total + 1
      if(red .eq. sc1_max) iblk0 = iblkmax
      ired(4) = iblk0
      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,
     &  MPI_MIN, MPI_COMM_WORLD, ierr)
      if(bl2pe(iblk0).eq.pe)then
        ired(1) = imaxsc1
        ired(2) = jmaxsc1
        ired(3) = kmaxsc1
        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,
     &    1, MPI_COMM_WORLD, sreq, ierr)
      endif
      if(pe.eq.0)then
        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &    1, MPI_COMM_WORLD, status, ierr)
        imaxsc1 = ired(1)
        jmaxsc1 = ired(2)
        kmaxsc1 = ired(3)
        iblkmax = ired(4)
      endif
      if(bl2pe(iblk0).eq.pe)then
        call MPI_WAIT(sreq,status,ierr)
      endif

      red = sc1_min
      call MPI_ALLREDUCE(red, sc1_min, 1, REAL_TYPE_MPI,
     &  MPI_MIN, MPI_COMM_WORLD, ierr)

      iblk0 = numbl_total + 1
      if(red .eq. sc1_min) iblk0 = iblkmin
      ired(4) = iblk0
      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,
     &  MPI_MIN, MPI_COMM_WORLD, ierr)
      if(bl2pe(iblk0).eq.pe)then
        ired(1) = iminsc1
        ired(2) = jminsc1
        ired(3) = kminsc1
        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,
     &    1, MPI_COMM_WORLD, sreq, ierr)
      endif
      if(pe.eq.0)then
        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),
     &    1, MPI_COMM_WORLD, status, ierr)
        iminsc1 = ired(1)
        jminsc1 = ired(2)
        kminsc1 = ired(3)
        iblkmin = ired(4)
      endif
      if(bl2pe(iblk0).eq.pe)then
        call MPI_WAIT(sreq,status,ierr)
      endif


      red = rho_max

      call MPI_ALLREDUCE(red, rho_max, 1, REAL_TYPE_MPI,

     &  MPI_MAX, MPI_COMM_WORLD, ierr)



      iblk0 = numbl_total + 1

      if(red .eq. rho_max) iblk0 = iblkmaxrho

      ired(4) = iblk0

      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,

     &  MPI_MIN, MPI_COMM_WORLD, ierr)

      if(bl2pe(iblk0).eq.pe)then

        ired(1) = imaxrho

        ired(2) = jmaxrho

        ired(3) = kmaxrho

        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,

     &    1, MPI_COMM_WORLD, sreq, ierr)

      endif

      if(pe.eq.0)then

        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),

     &    1, MPI_COMM_WORLD, status, ierr)

        imaxrho = ired(1)

        jmaxrho = ired(2)

        kmaxrho = ired(3)

        iblkmaxrho = ired(4)

      endif

      if(bl2pe(iblk0).eq.pe)then

        call MPI_WAIT(sreq,status,ierr)

      endif



      red = rho_min

      call MPI_ALLREDUCE(red, rho_min, 1, REAL_TYPE_MPI,

     &  MPI_MIN, MPI_COMM_WORLD, ierr)



      iblk0 = numbl_total + 1

      if(red .eq. rho_min) iblk0 = iblkminrho

      ired(4) = iblk0

      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,

     &  MPI_MIN, MPI_COMM_WORLD, ierr)

      if(bl2pe(iblk0).eq.pe)then

        ired(1) = iminrho

        ired(2) = jminrho

        ired(3) = kminrho

        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,

     &    1, MPI_COMM_WORLD, sreq, ierr)

      endif

      if(pe.eq.0)then

        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),

     &    1, MPI_COMM_WORLD, status, ierr)

        iminrho = ired(1)

        jminrho = ired(2)

        kminrho = ired(3)

        iblkminrho = ired(4)

      endif

      if(bl2pe(iblk0).eq.pe)then

        call MPI_WAIT(sreq,status,ierr)

      endif



c      red = T_max

c      call MPI_ALLREDUCE(red, T_max, 1, REAL_TYPE_MPI,

c     &  MPI_MAX, MPI_COMM_WORLD, ierr)



c      iblk0 = numbl_total + 1

c      if(red .eq. T_max) iblk0 = iblkmaxT

c      ired(4) = iblk0

c      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,

c     &  MPI_MIN, MPI_COMM_WORLD, ierr)

c      if(bl2pe(iblk0).eq.pe)then

c        ired(1) = imaxT

c        ired(2) = jmaxT

c        ired(3) = kmaxT

c        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,

c     &    1, MPI_COMM_WORLD, sreq, ierr)

c      endif

c      if(pe.eq.0)then

c        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),

c     &    1, MPI_COMM_WORLD, status, ierr)

c        imaxT = ired(1)

c        jmaxT = ired(2)

c        kmaxT = ired(3)

c        iblkmaxT = ired(4)

c      endif

c      if(bl2pe(iblk0).eq.pe)then

c        call MPI_WAIT(sreq,status,ierr)

c      endif



c      red = T_min

c      call MPI_ALLREDUCE(red, T_min, 1, REAL_TYPE_MPI,

c     &  MPI_MIN, MPI_COMM_WORLD, ierr)



c      iblk0 = numbl_total + 1

c      if(red .eq. T_min) iblk0 = iblkminT

c      ired(4) = iblk0

c      call MPI_ALLREDUCE(ired(4), iblk0, 1, INTEGER_TYPE_MPI,

c     &  MPI_MIN, MPI_COMM_WORLD, ierr)

c      if(bl2pe(iblk0).eq.pe)then

c        ired(1) = iminT

c        ired(2) = jminT

c        ired(3) = kminT

c        call MPI_ISEND(ired, 4, INTEGER_TYPE_MPI, 0,

c     &    1, MPI_COMM_WORLD, sreq, ierr)

c      endif

c      if(pe.eq.0)then

c        call MPI_RECV(ired, 4, INTEGER_TYPE_MPI, bl2pe(iblk0),

c     &    1, MPI_COMM_WORLD, status, ierr)

c        iminT = ired(1)

c        jminT = ired(2)

c        kminT = ired(3)

c        iblkminT = ired(4)

c      endif

c      if(bl2pe(iblk0).eq.pe)then

c        call MPI_WAIT(sreq,status,ierr)

c      endif


      red = u_max
      call MPI_REDUCE(red, u_max, 1, REAL_TYPE_MPI,
     &  MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      red = v_max
      call MPI_REDUCE(red, v_max, 1, REAL_TYPE_MPI,
     &  MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      red = w_max
      call MPI_REDUCE(red, w_max, 1, REAL_TYPE_MPI,
     &  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      red = u_min
      call MPI_REDUCE(red, u_min, 1, REAL_TYPE_MPI,
     &  MPI_MIN, 0, MPI_COMM_WORLD, ierr)
      red = v_min
      call MPI_REDUCE(red, v_min, 1, REAL_TYPE_MPI,
     &  MPI_MIN, 0, MPI_COMM_WORLD, ierr)
      red = w_min
      call MPI_REDUCE(red, w_min, 1, REAL_TYPE_MPI,
     &  MPI_MIN, 0, MPI_COMM_WORLD, ierr)

      if (master) then
        write(*,'(''CFL number: '',g14.8, '' at '',3i5,
     &    ''   in block '',i5)')
     &    cflmax,imxcfl,jmxcfl,kmxcfl,bmxcfl
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SAVE MAXDFS FOR ADJUSTING DELTA_T

      maxDFS = dfsmax

!        write(*,'(''Viscous stability number: '',g14.8, '' at '',3i5,
!     &    ''   in block '',i5)')
!     &    dfsmax,imxdfs,jmxdfs,kmxdfs,bmxdfs
        write(*,'(''Highest velocity: '',2(g14.8,'',''),g14.8)')
     &    u_max, v_max, w_max
        write(*,'(''Lowest velocity: '',2(g14.8,'',''),g14.8)')
     &    u_min, v_min, w_min
!        write(*,'(''The maximum scalar: '',g14.8, '' at '',3i5,
!     &    ''   in block '',i5)')
!     &    sc1_max,imaxsc1,jmaxsc1,kmaxsc1,iblkmax
!        write(*,'(''The minimum scalar: '',g14.8, '' at '',3i5,
!     &    ''   in block '',i5)')
!     &    sc1_min,iminsc1,jminsc1,kminsc1,iblkmin


        write(*,'(''Highest density: '',g14.8, '' recorded at '',3i5,

     &    ''   in block '',i5)')
     &    rho_max,imaxrho,jmaxrho,kmaxrho,iblkmaxrho
        write(*,'(''Lowest density: '',g14.8, '' recorded at '',3i5,

     &    ''   in block '',i5)')
     &    rho_min,iminrho,jminrho,kminrho,iblkminrho
c        write(*,'(''The maximum Temperature: '',g14.8, '' at '',3i5,

c     &    ''   in block '',i5)')
c     &    T_max,imaxT,jmaxT,kmaxT,iblkmaxT
c        write(*,'(''The minimum Temperature: '',g14.8, '' at '',3i5,

c     &    ''   in block '',i5)')
c     &    T_min,iminT,jminT,kminT,iblkminT

      end if

      return
      end
