!----------------------------------------------
! Step handles all of the simulation procedures
!----------------------------------------------

      subroutine step(ia,ra,u,v,w,p,sc1,
     &  sc1var,sc1dis,sc1varm,sc1dism,
     &  u_int,v_int,w_int,
     &  hn_u,hn_v,hn_w,hn_sc1,
     &  aw,ae,as,an,al,ar,ap,
     &  s,rho,res,gam,gamsc1,gamlam,
     &  x,y,z,
     &  dip,diu,div,dix,rdip,rdiv,rdix,
     &  djp,dju,djv,djx,rdjp,rdju,rdjx,
     &  hijp,hjip,hkip,hkjp,hiju,hjiu,hkiu,hkju,hijv,hjiv,hkiv,hkjv,
     &  rp,ru,rv,
     &  didjp,rdidjp,didju,rdidju,djdiv,rdjdiv,
     &  dkp,dkw,
     &  tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &  tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &  tw_wfk,tw_efk,tw_sfk,tw_nfk,rhonp1,rhonm1,Drhonp1Dt,  
     &  Temperature,
     &  srhot,srhoux, s_update)

      implicit none


#include "mpif.h"
#include "dim.inc"
#include "comm.inc"
#include "les.inc"
#include "data_type.inc"
#include "cb.control.inc"
#include "cb.address.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc"
#include "bc_type.inc"


      integer ia(*)
      real    ra(*)
      real    u(*), v(*), w(*), p(*), sc1(*)
      real    sc1var(*), sc1dis(*)
      real    sc1varm(*), sc1dism(*)
      real    u_int(*), v_int(*), w_int(*)
      real    hn_u(*), hn_v(*), hn_w(*), hn_sc1(*)
      real    ap(*), aw(*), ae(*), as(*), an(*), ar(*), al(*)
      real    s(*), rho(*), res(*), gam(*), gamsc1(*), gamlam(*)
      real    x(*), y(*), z(*)
      real    dip(*), diu(*), div(*), dix(*), rdip(*), rdiv(*), rdix(*)
      real    djp(*), dju(*), djv(*), djx(*), rdjp(*), rdju(*), rdjx(*)
      real    hijp(*), hjip(*), hkip(*), hkjp(*)
      real    hiju(*), hjiu(*), hkiu(*), hkju(*)
      real    hijv(*), hjiv(*), hkiv(*), hkjv(*)
      real    rp(*), ru(*), rv(*)
      real    didjp(*), rdidjp(*)
      real    didju(*), rdidju(*)
      real    djdiv(*), rdjdiv(*)
      real    dkp(*), dkw(*)
      real    tw_sfi(*), tw_nfi(*), tw_lfi(*), tw_rfi(*),
     &        tw_wfj(*), tw_efj(*), tw_lfj(*), tw_rfj(*),
     &        tw_wfk(*), tw_efk(*), tw_sfk(*), tw_nfk(*)
      real    rhonp1(*), rhonm1(*), Drhonp1Dt(*)
      real    Temperature(*)
      real    srhot(*), srhoux(*), s_update(*)


! local parameters

      integer     iseed, ik, io, jo, ko
      integer     ijk, ijks, ijk1, ijk2, k, j, i, ibl, ijkpai
      integer     imax, jmax, kmax,ijkn
      integer     pnt, ierr, ivar1, mpierr
      real        rdum(6)
      integer     idum(6)
      character*4 nmtd
      data        nmtd/'CENT'/

      real*4      etime, time_used(2)
      real        t_loop_start, t_loop_end, t_total
      real        uout

      integer itime

      logical     available

      real        rnorm2,fnorm2,enorm2,vnorm2
      ! Change these when the grid is changed.
      real        uprof(1:258), vprof(1:258),wprof(1:258)
      real        uuprof(1:130), vvprof(1:130),wwprof(1:130)
      real        umean(1:130),vmean(1:130),wmean(1:130)
      real        uumean(1:130),vvmean(1:130),wwmean(1:130)
       INTEGER, DIMENSION(:), ALLOCATABLE :: seed
       integer :: numb


! for debugging
      integer var_rotate
      common /debug1/var_rotate

      iseed = -30


! Set random number seed
      CALL RANDOM_SEED(SIZE=numb)
      ALLOCATE(seed(numb))

      OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
      READ(89) seed
      CLOSE(89)
      CALL RANDOM_SEED(put=seed)


! Set interfaces between blocks
      call set_interfaces

! Calculate metric quantities (dx, dy, dz and so on)
! Only done one as the grid is static
!       write(*,*) pe, 'starting comp metric'
!             call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call comp_metric(x,y,z,
     &     dip,diu,div,dix,rdip,rdiv,rdix,
     &     djp,dju,djv,djx,rdjp,rdju,rdjx,
     &     hijp,hjip,hkip,hkjp,
     &     hiju,hjiu,hkiu,hkju,
     &     hijv,hjiv,hkiv,hkjv,
     &     rp,ru,rv,
     &     didjp,rdidjp,didju,rdidju,djdiv,rdjdiv,
     &     dkp,dkw,
     &     ra(jw3a),ia(iwk1a))
!             call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!Open files for data output (flowviz, stats)
!One file per block      
      call open_files


!Initialise density arrays (fairly redundant but should be 
!tested before it can be removed)
      if(New_Run) then
         call init_array(rho,rholo,1,1,nummgl)
         call init_array(rhonp1,rholo,1,1,nummgl)
         call init_array(rhonm1,rholo,1,1,nummgl)
      endif

! Read in profiles for inflow condition
       if (master) write(*,*) 'reading inlet profiles'   
        call read_profile(uprof,vprof,wprof,uuprof,vvprof,wwprof)
 
! If this is a new simulation, initialise the flow field at t=0
      if(New_Run) then
       call initia(u,v,w,sc1,p,rho,rhonm1,rhonp1,y,uprof,uuprof,
     &              vprof,vvprof,wwprof) 
      endif
      
! read stats back in (need addback to be equal to 1 in input.dat file
! wam 13/09/13      
      if(New_Run .eq. .false. .and. addback .eq. 1) then
      call read_stats(ra(ju_m),ra(jv_m),ra(jw_m),
     &      ra(jruv),ra(jrvw),ra(jrwu),ra(jruu),ra(jrvv),ra(jrww),
     &      ra(jrrh),ra(jf_m),ra(jrff),
     &      ra(jTemperature), ra(jT_m),ra(jp_m),ra(jrho_m),ra(jsc1varm))
      endif      

!Time advancement is performed using the Adams-Bashforth scheme

      t_total = 0.
      iststep = 1
      phytime = 0.0

! Start of main time-stepping loop here
      do  istep = step_no, last
         if(master) then  
         write(*,*) '############ Time Step= ', istep,' ############'
         endif

!Add one to the iststep counter
!This is used in the calculation of the wall shear stress
!This could be brought into line with the way that the statistics
!are computed in statistics.F if needed.

      iststep = iststep + 1

! Set the timer to compute the time taken per time step
        t_loop_start = MPI_WTIME()


       if (master) then
         if (istep .eq. 1) then
          write(*,'(''Time step set to: '',i10)') dtime 
         endif
         

! Compute the relaxation factor for the density field
      alpha = exp(-1.0*dtime/t_relax)

      if (istep .eq. 1) then
      write(*,'(''Relaxation parameter alpha: '',g15.8)') alpha
	  endif
	  
       end if !master
        

! Send the time step and relaxation information to all processes
        call MPI_BCAST(dtime, 1, REAL_TYPE_MPI, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(phytime,1, REAL_TYPE_MPI,0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha,1, REAL_TYPE_MPI,0, MPI_COMM_WORLD, ierr) 

! Compute each solution variable in turn
! This loop solves for u,v,w,p,sc1 in turn
! the if statements in the main loop select the
! correct variable to solve.  This is rather inelegant
! and can be improved.

        do ivar = 1,numsolve

          gotu = ivar .eq. SV_U
          gotv = ivar .eq. SV_V
          gotw = ivar .eq. SV_W
          gotp = ivar .eq. SV_P
          gotsc1 = ivar .eq. SV_SC1

          call set_loop_bound(ivar,1,nummgl,1,numbl,2,idum)

! Set the inlet condition for the simulation
          if (gotu) then
            call fixed_inlet(hn_u,hn_v,hn_w,u_int,v_int,w_int,
     &        u,v,w,p,sc1,rho,gam,iseed,
     &        y,uprof,vprof,wprof,uuprof,vvprof,wwprof,umean,
     &        vmean,wmean,uumean,vvmean,wwmean,istep)
          end if
          
! Calculate the eddy viscosity
! The van-driest damping function is also computed here

          if (gotu) then
            call set_virt_bc(SV_U,1,u,ra(ju_m),rho,gam,
     &           ra(jw3a),ia(iwk1a),0.,
     &           tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &           tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &           tw_wfk,tw_efk,tw_sfk,tw_nfk,
     &           dip,diu,div,dix,djp,dju,djv,djx,ru,rv,dkp,dkw )

            call set_virt_bc(SV_V,1,v,ra(jv_m),rho,gam,
     &           ra(jw3a),ia(iwk1a),0.,
     &           tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &           tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &           tw_wfk,tw_efk,tw_sfk,tw_nfk,
     &           dip,diu,div,dix,djp,dju,djv,djx,ru,rv,dkp,dkw )

            call set_virt_bc(SV_W,1,w,ra(jw_m),rho,gam,
     &           ra(jw3a),ia(iwk1a),0.,
     &           tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &           tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &           tw_wfk,tw_efk,tw_sfk,tw_nfk,
     &           dip,diu,div,dix,djp,dju,djv,djx,ru,rv,dkp,dkw )
          end if
          if (gotu) then
            call subgrid(gam,gamsc1,u,v,w,sc1,rho,sc1var,sc1dis,
     &           ra(jw3a),dip,diu,djp,djv,dkp,dkw,
     &           hijp,hjip,hkip,hkjp,rp,ra(jw1a),ia(iwk1a))
!          write(*,*) pe,'eddy viscosity set'
          end if


! Solve the passive scalar field

          if (gotu) then
            call set_gamsc1_bc(sc1,gamsc1,gamlam,diu,djv,dkp,rp)

            call set_sc1(u,v,w,p,sc1,rho,gamsc1,gamlam,hn_sc1,
     &                   dip,diu,rdiv,djp,djv,rdju,
     &                   didjp,rdidjp,
     &                   hkip,hkjp,hkiu,hkjv,
     &                   rp,dkp,dkw,
     &                   ra(jw1a),ia(iwk1a))

            call set_sc1_bc(sc1)



           call  set_sc1var(u,v,w,sc1,rho,sc1var,sc1dis,
     &                      dip,djp,dkp,dkw,rp,ra(jw1a),ia(iwk1a))
     

! Update the density field from the scalar information
! Again this is rather inelegant

           call search_lib(
     &                    sc1, sc1var, rhonp1, Temperature,
     &                    ra(jw1a),ia(iwk1a),istep,
     &                    alpha, rho,rholo,rhohi)


! Set boundary conditions for the updated scalar field
            call set_sc1_bc(sc1)
            call set_sc1_bc(rhonp1)
            call set_sc1_bc(Temperature)


          endif


!Compute the intermediate velocities

          if (gotu) then

            call set_uv_int(ivar,u,v,w,p,rho,rhonp1,
     &           gam,u_int,hn_u,
     &           dip,diu,dix,rdix,djp,dju,djv,djx,rdjp,rdju,
     &           didjp,didju,rdidju,
     &           hijp,hiju,hjiu,hkiu,hkju,
     &           rp,ru,dkp,dkw,
     &           tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &           ra(jw1a),ia(iwk1a))


          else if (gotv) then

            call set_uv_int(ivar,v,u,w,p,rho,rhonp1,
     &           gam,v_int,hn_v,
     &           djp,djv,djx,rdjx,dip,div,diu,dix,rdip,rdiv,
     &           didjp,djdiv,rdjdiv,
     &           hjip,hjiv,hijv,hkjv,hkiv,
     &           rp,rv,dkp,dkw,
     &           tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &           ra(jw1a),ia(iwk1a))
        
     
          else if (gotw) then

            call set_w_int(u,v,w,p,rho,rhonp1,
     &           gam,w_int,hn_w,
     &           dip,diu,rdiv,djp,djv,rdju,
     &           didjp,rdidjp,
     &           hkip,hkjp,hkiu,hkjv,
     &           rp,ru,dkp,dkw,
     &           tw_sfk,tw_nfk,tw_wfk,tw_efk,
     &           ra(jw1a),ia(iwk1a))
          end if


!Apply the convective outflow condition
         
          call conv_exitbc(1,ivar,rho,rhonp1,u,v,w,
     &         u_int,v_int,w_int,hn_u,hn_v,hn_w,
     &         dip,diu,dix,djp,djv,djx,
     &         rdiv,rdju,didjp,
     &         rp,ru,rv,dkp,dkw,
     &         rdidjp)
!         write(*,*) 'outflow set'  


          if (gotp) then

!Compute the coefficients for the pressure equation
!Coefficients are time invariant, so only needs to be done on the
!first time step

            if (istep .eq. 1) then
              call coeff_p(diu,djv,rdiv,rdju,didjp,rp,dkp,dkw,
     &             ap,aw,ae,as,an,al,ar)
            end if

! Obtain the source term for the pressure solver

            call sour_p(diu,djv,rdiv,rdju,didjp,dkp,rdidjp,dkw,
     &         u_int,v_int,w_int,rho,rhonp1,Drhonp1Dt,s,srhot,srhoux,
     &         ra(jw1a),ra(jw2a))

! Solve for pressure

            call solvers(ap,aw,ae,as,an,al,ar,p,s,res,ia,ra)

! Extrapolate the pressure field to the boundaries

            call extrap_p(p,ra(jw1a),ia(iwk1a))


! Subtract a reference pressure from the pressure field
            call p_ref_subtract(p)

!This is the centreline boundary condition for a polar grid
!Should be moved to slipbc.F
!            do ibl = 1, numbl
!
!               if (bl_bc(3,ibl) .eq. BCT_CENTRE_LINE) then

!                  imax = bl_nijkp1(1,ibl)
!                  jmax = bl_nijkp1(2,ibl)
!                  kmax = bl_nijkp1(3,ibl)

!                  do i = 1, imax
!                     do k = 1, kmax
!                        ijk = i + imax*(2-1) + imax*jmax*(k-1) +
!     &                        sg_off(ibl)
!                        ijks = ijk - imax
!                          if ( k .le. kmax/2) then
!                            ijkpai = ijk + imax*jmax* (kmax/2)
!                          else
!                            ijkpai = ijk - imax*jmax* (kmax/2)
!                        endif

!                        u(ijks) = u(ijk)
!                        u_int(ijks) = u_int(ijk)

!                        w(ijks) = -1.0 * w(ijk)
!                        w_int(ijks) = -1.0 * w_int(ijk)

!                        v(ijks) = 0.5*(v(ijk) + v(ijkpai))
!                        v_int(ijks) = 0.5*(v_int(ijk) + v_int(ijkpai))
!                     enddo
!                  enddo
!               endif
!            enddo


!Use the new pressure field to compute the actual velocity field 
            call uvw_adv(SV_U,diu,rp,dkw,u,u_int,rhonp1,p,
     &           ra(jw1a),ia(iwk1a),istep)

            call uvw_adv(SV_V,djv,rp,dkw,v,v_int,rhonp1,p,
     &           ra(jw1a),ia(iwk1a),istep)

            call uvw_adv(SV_W,diu,rp,dkw,w,w_int,rhonp1,p,
     &           ra(jw1a),ia(iwk1a),istep)

!Apply slip boundary conditions  (if specified)
      call slipbc(x,y,z,dip,diu,djp,djv,dkp,
     &             u,v,v_int,w,sc1,rho,p,istep)


!Apply the outflow condition to the actual velocity field
            do ivar1 = SV_U,SV_W

              call conv_exitbc(2,ivar1,rho,rhonp1,u,v,w,
     &             u_int,v_int,w_int,hn_u,hn_v,hn_w,
     &             dip,diu,dix,djp,djv,djx,
     &             rdiv,rdju,didjp,
     &             rp,ru,rv,dkp,dkw,
     &             rdidjp)  
            end do


!Compute useful parameters such as CFL and DFS number
            call cflest(u,v,w,sc1,rhonp1,Temperature,
     &                  dip,djp,dkw,rp,gam)


!Swap rhonp1 -> rho -> rhonm1 for next timestep
            call swap_rhos(rho,rhonp1,rhonm1)

          end if  !gotp

       enddo !ivar
       
!Accumulate statistical information
!These pointers are 3-d fields, whereas the code performs
!spanwise averaging, hence the outputted fields are 2-D.
!The pointers should be modified to save memory.

!Statistics are only accumlated after 'statstart'
!i.e. statstart should be the time step where the flow
!has reached a statistically stationary state

      if(istep .ge. statstart) then
       if(mod(istep, stat_freq) .eq. 0) then

        call statistics(u,v,w,p,rho,ra(ju_m),ra(jv_m),ra(jw_m),
     &      ra(jp_m), ra(jrho_m),
     &      ra(jruv),ra(jrvw),ra(jrwu),ra(jruu),ra(jrvv),ra(jrww),
     &      ra(jrrh),gam,gamsc1,
     &      ra(jsc1),ra(jf_m),ra(jrff),
     &      ra(jsc1var),ra(jsc1dis),ra(jsc1varm),ra(jsc1dism),
     &      ra(jTemperature), ra(jT_m),
     &      dip,djp,rp,dkp,dkw,
     &      x,y,z,istep)

        endif
      endif

!Compute wall time taken for the time step
        t_loop_end = MPI_WTIME()

        rdum(1) = t_loop_end-t_loop_start
        call MPI_REDUCE(rdum(1), rdum(2), 1, REAL_TYPE_MPI, MPI_MAX,
     &    0, MPI_COMM_WORLD, mpierr)

        if(master)then
          write(*,'(''Total time used for this step was'',g15.8,
     &  ''seconds''  )')
     &    rdum(2)
        end if
        t_total = t_total + t_loop_end-t_loop_start



!Output data to file

! Write out flow visualisation data
       if(istep .ge. fieldstart) then
        if(mod(istep,Fiel_freq).eq.0)then
       if(master) write(*,*) 'Writing out flow visualisation data'
      call field_out(x,y,z,dip,diu,djp,djv,dkp,
     &             u,v,w,sc1,rho,p,gam,istep)
        endif
        endif ! statstart

        

!Write out 3-D field data dumps
        if (mod(istep,dump_freq).eq.0 .or. istep.eq.last)then
          do ivar = 1,numvar
            if (ifile(ivar) .ne. 0) then
              pnt = pntvar(ivar)
              call file_3(ndn,namvar(ivar),wformat,ra(pnt),1,1, ierr)
            end if
           end do
        end if  !dump_freq
          

      if (mod(istep,dump_freq).eq.0 .or. istep.eq.last)then
          if (writewall) then
            call file_twall(tw_sfi,tw_nfi,tw_lfi,tw_rfi,
     &           tw_wfj,tw_efj,tw_lfj,tw_rfj,
     &           tw_wfk,tw_efk,tw_sfk,tw_nfk)
          end if
      endif

!Write out statistical information
      if(istep .ge. statstart) then
      if (mod(istep,dump_freq).eq.0 .or. istep.eq.last)then
           call dump_stats(x,y,z,ra(ju_m),ra(jv_m),ra(jw_m),
     &      ra(jruv),ra(jrvw),ra(jrwu),ra(jruu),ra(jrvv),ra(jrww),
     &      ra(jrrh),ra(jf_m),ra(jrff),
     &      ra(jTemperature), ra(jT_m),ra(jp_m),ra(jrho_m),ra(jsc1varm))
         end if
       endif ! istep
 
!Write out instantaneous data for the computation of spectra                
!       if( istep .ge. statstart) then
!       if (mod(istep,spec_freq).eq.0 .or. istep.eq.last)then
!         call spectra(x,y,z,dip,diu,djp,djv,dkp,
!     &             u,v,w,sc1,rho)
!       endif
!       endif

       enddo  !istep - end of main time step loop
       
      rdum(1) = t_total
      call MPI_ALLREDUCE(rdum, t_total, 1, REAL_TYPE_MPI, MPI_MAX,
     &  MPI_COMM_WORLD, mpierr)

      if(master)then
        write(*,'(''Calculation completed in'',g15.8,'' seconds'')')
     &  t_total
        close(87)
      endif

      if(master) then
      write(*,*) 'Final time step = ', istep
      endif

      return
      end subroutine
