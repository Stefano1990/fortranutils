!--------------------------------
! LESter 
!--------------------------------

      program LESter

      implicit none


#include "mpif.h"
#include "dim.inc"
#include "comm.inc"
#include "data_type.inc"


      integer i, ierr, iblk
      integer pvar
      logical die, available
      character*32 fsetup

      logical mpiflag
      integer iitemp, j, status(MPI_STATUS_SIZE)
      character*513 errstring

#include "cb.big.inc"
#include "cb.control.inc"
#include "cb.gridp.inc"
#include "cb.lgrid.inc"
#include "cb.address.inc"
#include "cb.comm.inc"
 
      data   die    /.false./

      integer var_rotate
      common /debug1/var_rotate

! Initialise large arrays
      if(master) then
      PRINT*, ' starting LESter'
      endif
      do i=1,maxra
        ra(i) = 0.0
      end do
      do i=1,maxia
        ia(i) = 0
      end do

! Set unit numbers for reading data

      ndn    = 90
      nblset = 92

! Set up parallelisation

      call setupcomm()


! Set number of variables and number of solution variables

      numvar   = 47
      numsolve = 5

! Read in the size of each computational block

      call read_bl_size(nblset)


!Define variable names (for outputting files)
!Note that not all files are outputted when dumping data
!only u,v,w,p,hn_u,hn_v,hn_w,sc1,rho and gam are needed
!to restart a calculation


      namvar(1 ) = 'u'
      namvar(2 ) = 'v'
      namvar(3 ) = 'w'
      namvar(4 ) = 'p'
      namvar(5 ) = 'sc1'
      namvar(6 ) = 'sc1var'
      namvar(7 ) = 'sc1dis'
      namvar(8 ) = 'u_int'
      namvar(9 ) = 'v_int'
      namvar(10 ) = 'w_int'
      namvar(11 ) = 'hn_u'
      namvar(12 ) = 'hn_v'
      namvar(13) = 'hn_w'
      namvar(14) = 'hn_sc1'

      namvar(15) = 'u_m'
      namvar(16) = 'v_m'
      namvar(17) = 'w_m'
      namvar(18) = 'f_m'

      namvar(19) = 'ruu'
      namvar(20) = 'rvv'
      namvar(21) = 'rww'
      namvar(22) = 'ruv'
      namvar(23) = 'rvw'
      namvar(24) = 'rwu'
      namvar(25) = 'rff'
      namvar(26) = 'rrh'
      namvar(27) = 'p_m'
      namvar(28) = 'rho_m' 
      namvar(29) = 'suu' 
      namvar(30) = 'svv' 
      namvar(31) = 'sww' 
      namvar(32) = 'suv' 
      namvar(33) = 'svw' 
      namvar(34) = 'swu' 
      namvar(35) = 'sc1varm'
      namvar(36) = 'sc1dism'
      namvar(37) = 'rho'
      namvar(38) = 'rhonp1'
      namvar(39) = 'rhonm1'
      namvar(40) = 'Drhonp1Dt'
      namvar(41) = 'Temperature'
      namvar(42) = 'gam'
      namvar(43) = 's'
      namvar(44) = 'srhot'
      namvar(45) = 'srhoux'
      namvar(46) = 's_update'
      namvar(47) = 'T_m' 


! Set default parameters
      call defaults
!These parameters can be overwritten if read in from input.dat

! Read simulation input parameters
      call dataread
      call read_bl_info(ierr)

!Send block information to relevant processes
      if(master)then
        do iblk=1,numbl_total
          if (pe.ne.bl2pe(iblk)) then
            call MPI_SEND(bl_nijkp1(1,iblk), 3, INTEGER_TYPE_MPI,
     &      bl2pe(iblk), 1, MPI_COMM_WORLD, ierr)
            bl_nijkp1_master(1,iblk) = bl_nijkp1(1,iblk)
            bl_nijkp1_master(2,iblk) = bl_nijkp1(2,iblk)
            bl_nijkp1_master(3,iblk) = bl_nijkp1(3,iblk)
          end if
        end do
        iitemp = 1
        do iblk=1,numbl_total
          if(pe.eq.bl2pe(iblk)) then
            bl_nijkp1_master(1,iblk) = bl_nijkp1(1,iblk)
            bl_nijkp1_master(2,iblk) = bl_nijkp1(2,iblk)
            bl_nijkp1_master(3,iblk) = bl_nijkp1(3,iblk)
            bl_nijkp1(1,iitemp) = bl_nijkp1(1,iblk)
            bl_nijkp1(2,iitemp) = bl_nijkp1(2,iblk)
            bl_nijkp1(3,iitemp) = bl_nijkp1(3,iblk)
            iitemp = iitemp + 1
          endif
        end do
      else
        do iblk=1,numbl
          call MPI_RECV(bl_nijkp1(1,iblk),3, INTEGER_TYPE_MPI, 0, 1,
     &    MPI_COMM_WORLD, status, ierr)
        end do
      endif

! Set the offsets
      call set_offset()

! Set some constants

      call set_const()

! Create array pointers

      call set_pointer()

      pntvar(1 ) = ju
      pntvar(2 ) = jv
      pntvar(3 ) = jw
      pntvar(4 ) = jp
      pntvar(5 ) = jsc1
      pntvar(6 ) = jsc1var
      pntvar(7 ) = jsc1dis
      pntvar(8 ) = ju_int
      pntvar(9 ) = jv_int
      pntvar(10 ) = jw_int
      pntvar(11 ) = jhn_u
      pntvar(12 ) = jhn_v
      pntvar(13) = jhn_w
      pntvar(14) = jhn_sc1

      pntvar(15) = ju_m
      pntvar(16) = jv_m
      pntvar(17) = jw_m
      pntvar(18) = jf_m

      pntvar(19) = jruu
      pntvar(20) = jrvv
      pntvar(21) = jrww
      pntvar(22) = jruv
      pntvar(23) = jrvw
      pntvar(24) = jrwu
      pntvar(25) = jrff
      pntvar(26) = jrrh
      pntvar(27) = jp_m
      pntvar(28) = jrho_m
      pntvar(29) = jsuu 
      pntvar(30) = jsvv 
      pntvar(31) = jsww 
      pntvar(32) = jsuv 
      pntvar(33) = jsvw 
      pntvar(34) = jswu       

      pntvar(35) = jsc1varm
      pntvar(36) = jsc1dism
      pntvar(37) = jrho
      pntvar(38) = jrhonp1
      pntvar(39) = jrhonm1
      pntvar(40) = jDrhonp1Dt

      pntvar(41) = jTemperature
      pntvar(42) = jgam
      pntvar(43) = js
      pntvar(44) = jsrhot
      pntvar(45) = jsrhoux
      pntvar(46) = js_update
      pntvar(47) = jT_m


!Read in the grid (can be cartesian or polar)


      call read_xyz(ndn,'x.bin',ra(jx),ra(jw1a),ierr)
      call read_xyz(ndn,'y.bin',ra(jy),ra(jw1a),ierr)
      call read_xyz(ndn,'z.bin',ra(jz),ra(jw1a),ierr)

! Construct the coarse grids for pressure solver

      call mg_xyz('x',ra(jx),ra(jw1a))
      call mg_xyz('y',ra(jy),ra(jw1a))
      call mg_xyz('z',ra(jz),ra(jw1a))

! Read in stored simulation data - 3d fields
      if(New_Run .eq. .false.) then
      if (ierr.eq.0) then
        do ivar = 1,numvar
            if (ifile(ivar) .ne. 0) then
          pvar = pntvar(ivar)
          call read_3(ndn,namvar(ivar),rformat,ra(pvar),ra(jw1a),ierr)
            endif
        end do
      end if
      endif

! Read in stored simulation data - 2-d fields
! Read in the wall shear stress
      if(.not. New_Run .and. writewall) then
      call read_twall(ra(jtw_sfi),ra(jtw_nfi),ra(jtw_lfi),ra(jtw_rfi),
     &     ra(jtw_wfj),ra(jtw_efj),ra(jtw_lfj),ra(jtw_rfj),
     &     ra(jtw_wfk),ra(jtw_efk),ra(jtw_sfk),ra(jtw_nfk),
     &     ra(jw1a),ierr)
      endif

!Check for any errors during the set-up process

      if (ierr.ne.0.and.master) then
        write(*,'(''LESter failed.  Please check parameters.'')')
        stop
      end if

! Start simulation procedure

      call step(ia,ra,ra(ju),ra(jv),ra(jw),ra(jp),ra(jsc1),
     &     ra(jsc1var),ra(jsc1dis),ra(jsc1varm),ra(jsc1dism),
     &     ra(ju_int),ra(jv_int),ra(jw_int),
     &     ra(jhn_u),ra(jhn_v),ra(jhn_w),ra(jhn_sc1),
     &     ra(jaw),ra(jae),ra(jas),ra(jan),ra(jal),ra(jar),ra(jap),
     &     ra(js),ra(jrho),ra(jres),ra(jgam),ra(jgamsc1),ra(jgamlam),
     &     ra(jx),ra(jy),ra(jz),ra(jdip),ra(jdiu),ra(jdiv),ra(jdix),
     &     ra(jrdip),ra(jrdiv),ra(jrdix),
     &     ra(jdjp),ra(jdju),ra(jdjv),ra(jdjx),
     &     ra(jrdjp),ra(jrdju),ra(jrdjx),
     &     ra(jhijp),ra(jhjip),ra(jhkip),ra(jhkjp),
     &     ra(jhiju),ra(jhjiu),ra(jhkiu),ra(jhkju),
     &     ra(jhijv),ra(jhjiv),ra(jhkiv),ra(jhkjv),
     &     ra(jrp),ra(jru),ra(jrv),
     &     ra(jdidjp),ra(jrdidjp),ra(jdidju),ra(jrdidju),
     &     ra(jdjdiv),ra(jrdjdiv),ra(jdkp),ra(jdkw),
     &     ra(jtw_sfi),ra(jtw_nfi),ra(jtw_lfi),ra(jtw_rfi),
     &     ra(jtw_wfj),ra(jtw_efj),ra(jtw_lfj),ra(jtw_rfj),
     &     ra(jtw_wfk),ra(jtw_efk),ra(jtw_sfk),ra(jtw_nfk),
     &     ra(jrhonp1),ra(jrhonm1),ra(jDrhonp1Dt),ra(jTemperature),
     &     ra(jsrhot),ra(jsrhoux),ra(js_update))

#if ! defined MPI1 && ( ! defined SHMEM)
      call MPI_WIN_FREE(workwin, ierr)
#endif
      call MPI_FINALIZE(ierr)

      stop
      end
