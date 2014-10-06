C=======================================================================
C APPLY THE CONVECTIVE BOUNDARY CONDITION TO THE EXIT BOUNDARY
C=======================================================================

C NOTE:
C 1. FOR VELOCITY ONLY
C 2. ONLY TWO TYPES OF BOUNDARY CONDITIONS ARE CONSIDERED TO BE NECES-
C    SARY TO CALCULATE MASS FLUX
C 3. BOUNDARIES WITH B.C. CONDITION BCT_FIXED ARE TAKEN TO BE INLET
C    BOUNDARIES
C 4. TWO TASKS ARE IDENTIFIED HERE:
C    A. task = 1 FOR CALCULATING INTERMEDIATE VELOCITY
C    B. task = 2 FOR CALCULATING THE VELOCITY AT THE NEW TIME STEP n+1
C 5. IT IS ASSUMED THAT THE INITIAL VELOCITY FIELD GIVES NON-ZERO
C    FLUXES.
C 6. IN mflux_out, THE CONTRIBUTIONS FROM BOUNDARIES WITH OTHER B.C.
C    THAN THE CONVECTIVE B.C. HAVE NOT BEEN TAKEN INTO ACCOUNT; THESE
C    CONTRIBUTIONS SHOULD BE ADDED TO mflux_out.


      subroutine conv_exitbc(task,iva,rho,rhonp1,u,v,w,
     &  u_int,v_int,w_int,hu,hv,hw,
     &  dip,diu,dix,djp,djv,djx,
     &  rdiv,rdju,didjp,
     &  rp,ru,rv,dkp,dkw,
     &  rdidjp)  

      implicit none

C MPI
#include "mpif.h"

C PARAMETERS
#include "dim.inc"
#include "bc_type.inc"
#include "les.inc"
#include "comm.inc"
#include "data_type.inc"

C COMMON BLOCKS
#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc"


      integer task, iva
      real    rho(*)
      real    u(*), v(*), w(*)
      real    u_int(*), v_int(*), w_int(*)
      real    hu(*), hv(*), hw(*)
      real    dip(*), diu(*), dix(*), djp(*), djv(*), djx(*)
      real    rdiv(*), rdju(*), didjp(*)
      real    rp(*), ru(*), rv(*)
      real    dkp(*), dkw(*)
      real    rdidjp(*)  

      real    rhonp1(*)

C LOCAL

      integer p1, p2, p3
      integer nip1, njp1, nkp1
      integer io, jo, ko
      integer d1, iside, iblk, iplane, iva1
      integer ijkse(6), ps_ijk
      integer istr, iend, jstr, jend, kstr, kend
      integer mod_iside
      real    mflux, mflux_in, mflux_out
      real    area, uave, scale_factor, area_rho
      logical first_call, no_convective_bc, ltemp
      real    mchg_rate 
      integer ierr
      integer status(MPI_STATUS_SIZE), i
      real    red, red2(6),uavetemp
#ifdef DEBUG
      real    mflux_out2(maxbl), mflux_in2(maxbl)
#endif

      data    first_call/.true./

      save    mflux_in, mflux_out
      save    first_call, no_convective_bc

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C THE FIRST CALL OF THIS ROUTINE WILL DETERMINE WHETHER THERE ARE
C CONVECTIVE BOUNDARY CONDITIONS. IF THERE ARE NOT, IT RETURNS
C IMMEDIATELY IN THE FOLLOWING CALLS.
C-----------------------------------------------------------------------

      if (first_call) then
        first_call       = .false.
        no_convective_bc = .true.
        do iblk = 1,numbl
          do iside = 1,6
            if (bl_bc(iside,iblk) .eq. BCT_CONVECTIVE) then
              no_convective_bc = .false.
              goto 5
            end if
          end do
        end do
 5      continue
        ltemp = no_convective_bc
        call MPI_ALLREDUCE(ltemp, no_convective_bc, 1, MPI_LOGICAL,
     &  MPI_LAND, MPI_COMM_WORLD, ierr)
      else
        if (no_convective_bc) return
      end if

      if (no_convective_bc) return

C-----------------------------------------------------------------------
C THE OUTFLOW MASS FLUXES AT EXIT(S) (USING VELOCITY AT STEP n). THESE
C MASS FLUXES ARE USED TO CALCULATE THE AVERAGE VELOCITY AT EXITS
C-----------------------------------------------------------------------

      if (task.eq.1 .and. iva.eq.SV_U) then
        
        do iside = 1,6	  
            mfluxbl_out_ttl(iside) =0.0
            area_rho_exitbl_ttl(iside) =0.0
        enddo

        do iblk = 1,numbl
          do iside = 1,6

            if (bl_bc(iside,iblk) .eq. BCT_CONVECTIVE) then

              p1 = mg_off_1d_k(iblk,1) + 1
              p2 = mg_off_2d_ij(iblk,1) + 1
              p3 = mg_off(iblk,1) + 1

              nip1 = mg_nijkp1(1,iblk,1)
              njp1 = mg_nijkp1(2,iblk,1)
              nkp1 = mg_nijkp1(3,iblk,1)

              d1 = (iside+1)/2

              call set_loop_bound(d1,1,1,iblk,iblk,1,ijkse)

              mod_iside = mod(iside,2)
              iplane = (ijkse(2*d1-1)-1)*mod_iside +
     &                 (ijkse(2*d1  )+1)*(1-mod_iside)

              call mass_flux(iplane,nip1,njp1,nkp1,d1,d1,
     &             ijkse(1),ijkse(2),ijkse(3),ijkse(4),
     &             ijkse(5),ijkse(6),

C  !!! the bulk velocity in task 1 and 2 are all calculated using rhonp1
     &             area,area_rho,mflux,u(p3),v(p3),w(p3),rhonp1(p3),
     &             rdju(p2),rdiv(p2),didjp(p2),dkp(p1))

              area_exitbl(iside,iblk) = area
              mfluxbl_out(iside,iblk) = mflux
              area_rho_exitbl(iside,iblk) = area_rho

Cwk The following treatment is ONLY suitable for single exit conditon
              mfluxbl_out_ttl(iside) = mfluxbl_out_ttl(iside) + 
     &                          mfluxbl_out(iside,iblk)
              area_rho_exitbl_ttl(iside) = area_rho_exitbl_ttl(iside) + 
     &                              area_rho_exitbl(iside,iblk)
            end if
          end do
        end do

        do iside = 1,6        
           red2(iside) = mfluxbl_out_ttl(iside)
        enddo
        call MPI_ALLREDUCE(red2, mfluxbl_out_ttl, 6, 
     &  REAL_TYPE_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)

        do iside = 1,6        
           red2(iside) = area_rho_exitbl_ttl(iside)
        enddo
        call MPI_ALLREDUCE(red2, area_rho_exitbl_ttl, 6, 
     &  REAL_TYPE_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)

      end if

C-----------------------------------------------------------------------
C MIMIC THE ADAMS-BASHFORTH SCHEME TO OBTAIN THE INTERMEDIATE VELOCITY,
C I.E., USE FORMULA:
C                       n                 n-1
C  ^    n   3dt   n pd(u )   dt   n-1 pd(u   )
C  u = u  - ---*ue  ------ + --*ue    --------
C            2       pd x    2          pd x
C-----------------------------------------------------------------------

      if (iva .eq. SV_U) then
        mflux_in  = 0.0
        mflux_out = 0.0
      end if
#ifdef DEBUG
      do iblk=1,numbl
        mflux_in2(iblk) = 0.0
        mflux_out2(iblk) = 0.0
      end do
#endif

      if (task .eq. 1) then

        do iblk = 1,numbl
          do iside = 1,6

C - THE TOTAL INFLOW MASS FLUX (USING THE INTERMEDIATE VELOCITY) (THIS
C   IS USED IN SCALING)

            p1 = mg_off_1d_k(iblk,1) + 1
            p2 = mg_off_2d_ij(iblk,1) + 1
            p3 = mg_off(iblk,1) + 1

            nip1 = mg_nijkp1(1,iblk,1)
            njp1 = mg_nijkp1(2,iblk,1)
            nkp1 = mg_nijkp1(3,iblk,1)

            d1 = (iside+1)/2

            call set_loop_bound(d1,1,1,iblk,iblk,1,ijkse)
            mod_iside = mod(iside,2)


            iplane = (ijkse(2*d1-1)-1)*mod_iside +
     &               (ijkse(2*d1  )+1)*(1-mod_iside)

            if (bl_bc(iside,iblk) .eq. BCT_FIXED) then
              call mass_flux(iplane,nip1,njp1,nkp1,d1,iva,
     &             ijkse(1),ijkse(2),ijkse(3),ijkse(4),
     &             ijkse(5),ijkse(6),
     &             area,area_rho,
cwk     &             mflux,u_int(p3),v_int(p3),w_int(p3),rho(p3),
     &             mflux,u_int(p3),v_int(p3),w_int(p3),rhonp1(p3),
     &             rdju(p2),rdiv(p2),didjp(p2),dkp(p1))
              mflux_in = mflux_in + mflux
#ifdef DEBUG
              mflux_in2(iblk) = mflux_in2(iblk) + mflux
#endif
            end if

            if (bl_bc(iside,iblk) .eq. BCT_CONVECTIVE) then

              io = 1
              jo = nip1
              ko = nip1*njp1

C - INTERMEDIATE VELOCITY

              uave=mfluxbl_out_ttl(iside)/area_rho_exitbl_ttl(iside)        



!         write(*,1001)  mfluxbl_out_ttl(iside),
!     &                       area_rho_exitbl_ttl(iside), uave, iblk,iva
! 1001         format('*t1* mflux_ttl= ',e10.4,' arho_ttl= ',e10.4,
!     &        ' uave= ',e10.4, ' block= ', i4, ' iva= ', i4)

              call set_loop_bound(iva,1,1,iblk,iblk,1,ijkse)
              istr = ijkse(1)
              iend = ijkse(2)
              jstr = ijkse(3)
              jend = ijkse(4)
              kstr = ijkse(5)
              kend = ijkse(6)


              if (d1 .eq. SV_U) then

                if (iva .eq. SV_U) then
                  call conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
cwk     &                 u(p3),u_int(p3),rho(p3),hu(p3),
     &                 u(p3),u_int(p3),rho(p3),rhonp1(p3),hu(p3),
     &                 dip(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_V) then
                  call conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
cwk     &                 v(p3),v_int(p3),rho(p3),hv(p3),
     &                 v(p3),v_int(p3),rho(p3),rhonp1(p3),hv(p3),
     &                 dix(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_W) then
                  call conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
cwk     &                 w(p3),w_int(p3),rho(p3),hw(p3),
     &                 w(p3),w_int(p3),rho(p3),rhonp1(p3),hw(p3),
     &                 diu(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                end if

              else if (d1 .eq. SV_V) then

                if (iva .eq. SV_U) then
                  call conv_exitbc1_block(d1,iside,njp1,nip1,nkp1,
     &                 jstr,jend,istr,iend,kstr,kend,jo,io,ko,
cwk     &                 u(p3),u_int(p3),rho(p3),hu(p3),
     &                 u(p3),u_int(p3),rho(p3),rhonp1(p3),hu(p3),
     &                 djx(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_V) then
                  call conv_exitbc1_block(d1,iside,njp1,nip1,nkp1,
     &                 jstr,jend,istr,iend,kstr,kend,jo,io,ko,
cwk     &                 v(p3),v_int(p3),rho(p3),hv(p3),
     &                 v(p3),v_int(p3),rho(p3),rhonp1(p3),hv(p3),
     &                 djp(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_W) then
                  call conv_exitbc1_block(d1,iside,njp1,nip1,nkp1,
     &                 jstr,jend,istr,iend,kstr,kend,jo,io,ko,
cwk     &                 w(p3),w_int(p3),rho(p3),hw(p3),
     &                 w(p3),w_int(p3),rho(p3),rhonp1(p3),hw(p3),
     &                 djv(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                end if

              else

                if (iva .eq. SV_U) then
                  call conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
cwk     &                 u(p3),u_int(p3),rho(p3),hu(p3),
     &                 u(p3),u_int(p3),rho(p3),rhonp1(p3),hu(p3),
     &                 djx(p2),ru(p2),dkw(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_V) then
                  call conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
cwk     &                 v(p3),v_int(p3),rho(p3),hv(p3),
     &                 v(p3),v_int(p3),rho(p3),rhonp1(p3),hv(p3),
     &                 djx(p2),rv(p2),dkw(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_W) then
                  call conv_exitbc1_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
CWK     &                 w(p3),w_int(p3),rho(p3),hw(p3),
     &                 w(p3),w_int(p3),rho(p3),rhonp1(p3),hw(p3),
     &                 djx(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                end if

              end if

C - THE TOTAL OUTFLOW MASS FLUX (USING THE INTERMEDIATE VELOCITY) (THIS
C   IS USED IN SCALING)

              call mass_flux(iplane,nip1,njp1,nkp1,d1,iva,
     &             ijkse(1),ijkse(2),ijkse(3),ijkse(4),
     &             ijkse(5),ijkse(6),
     &             area,area_rho,
cwk     &             mflux,u_int(p3),v_int(p3),w_int(p3),rho(p3),
     &             mflux,u_int(p3),v_int(p3),w_int(p3),rhonp1(p3),
     &             rdju(p2),rdiv(p2),didjp(p2),dkp(p1))
              mflux_out = mflux_out + mflux
              mfluxbl_out2(iside,iblk) = mflux +
     &                                   mfluxbl_out2(iside,iblk)
#ifdef DEBUG
              mflux_out2(iblk) = mflux_out2(iblk) + mflux
#endif

            end if

          end do
        end do

C-----------------------------------------------------------------------
C                                                        n+1
C USE THE CONDITION DIRECTLY TO OBTAIN THE NEW VELOCITY u    FOR THE
C EXIT BOUNDARY, I.E.,
C    n+1    n          n+1    n+1
C   u(i) - u(i)     n u(i) - u(i-1)
C   ----------- + ue  ------------- = 0
C       dt                 dx
C----------------------------------------------------------------------

      else

        do iblk = 1,numbl
          do iside = 1,6

            p1 = mg_off_1d_k(iblk,1) + 1
            p2 = mg_off_2d_ij(iblk,1) + 1
            p3 = mg_off(iblk,1) + 1

            nip1 = mg_nijkp1(1,iblk,1)
            njp1 = mg_nijkp1(2,iblk,1)
            nkp1 = mg_nijkp1(3,iblk,1)

            d1 = (iside+1)/2

            call set_loop_bound(d1,1,1,iblk,iblk,1,ijkse)
            mod_iside = mod(iside,2)

            iplane = (ijkse(2*d1-1)-1)*mod_iside +
     &               (ijkse(2*d1  )+1)*(1-mod_iside)

C - THE TOTAL INFLOW MASS FLUX (USED IN SCALING)

            if (bl_bc(iside,iblk) .eq. BCT_FIXED) then
              call mass_flux(iplane,nip1,njp1,nkp1,d1,iva,
     &             ijkse(1),ijkse(2),ijkse(3),ijkse(4),
     &             ijkse(5),ijkse(6),
     &             area,area_rho,mflux,u(p3),v(p3),w(p3),rhonp1(p3),
     &             rdju(p2),rdiv(p2),didjp(p2),dkp(p1))
              mflux_in = mflux_in + mflux
#ifdef DEBUG
              mflux_in2(iblk) = mflux_in2(iblk) + mflux
#endif
            end if

            if (bl_bc(iside,iblk) .eq. BCT_CONVECTIVE) then

              io = 1
              jo = nip1
              ko = nip1*njp1

C - VELCOTIY VELOCITY AT THE NEW TIME STEP


C++ NOTE: uave is not a common value, thus has to be calculated here for task2
              uave=mfluxbl_out_ttl(iside)/area_rho_exitbl_ttl(iside)      


!     &                       area_rho_exitbl_ttl(iside), uave, iblk, iva
! 1002         format('*t2* mflux_ttl= ',e10.4,' arho_ttl= ',e10.4,
!     &        ' uave= ',e10.4, ' block= ', i4, ' iva=', i4)

              call set_loop_bound(iva,1,1,iblk,iblk,1,ijkse)
              istr = ijkse(1)
              iend = ijkse(2)
              jstr = ijkse(3)
              jend = ijkse(4)
              kstr = ijkse(5)
              kend = ijkse(6)

              if (d1 .eq. SV_U) then

                if (iva .eq. SV_U) then
                  call conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &                 u(p3),u_int(p3),
     &                 rhonp1(p3),dip(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_V) then
                  call conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &                 v(p3),v_int(p3),
     &                 rhonp1(p3),dix(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_W) then
                  call conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &                 w(p3),w_int(p3),
     &                 rhonp1(p3),diu(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                end if

              else if (d1 .eq. SV_V) then

                if (iva .eq. SV_U) then
                  call conv_exitbc2_block(d1,iside,njp1,nip1,nkp1,
     &                 jstr,jend,istr,iend,kstr,kend,jo,io,ko,
     &                 u(p3),u_int(p3),
     &                 rhonp1(p3),djx(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_V) then
                  call conv_exitbc2_block(d1,iside,njp1,nip1,nkp1,
     &                 jstr,jend,istr,iend,kstr,kend,jo,io,ko,
     &                 v(p3),v_int(p3),
     &                 rhonp1(p3),djp(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_W) then
                  call conv_exitbc2_block(d1,iside,njp1,nip1,nkp1,
     &                 jstr,jend,istr,iend,kstr,kend,jo,io,ko,
     &                 w(p3),w_int(p3),
     &                 rhonp1(p3),djv(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                end if

              else

                if (iva .eq. SV_U) then
                  call conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &                 u(p3),u_int(p3),
     &                 rhonp1(p3),djx(p2),ru(p2),dkw(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_V) then
                  call conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &                 v(p3),v_int(p3),
     &                 rhonp1(p3),djx(p2),rv(p2),dkw(p1),
     &                 dtime,uave)
                else if (iva .eq. SV_W) then
                  call conv_exitbc2_block(d1,iside,nip1,njp1,nkp1,
     &                 istr,iend,jstr,jend,kstr,kend,io,jo,ko,
     &                 w(p3),w_int(p3),
     &                 rhonp1(p3),djx(p2),rp(p2),dkp(p1),
     &                 dtime,uave)
                end if

              end if

C - THE TOTAL OUTFLOW MASS FLUX (USED IN SCALING)

              call mass_flux(iplane,nip1,njp1,nkp1,d1,iva,
     &             ijkse(1),ijkse(2),ijkse(3),ijkse(4),
     &             ijkse(5),ijkse(6),
cwk     &             area,area_rho,mflux,u(p3),v(p3),w(p3),rho(p3),
     &             area,area_rho,mflux,u(p3),v(p3),w(p3),rhonp1(p3),
     &             rdju(p2),rdiv(p2),didjp(p2),dkp(p1))
              mflux_out = mflux_out + mflux
              mfluxbl_out2(iside,iblk) = mflux +
     &                                   mfluxbl_out2(iside,iblk)
#ifdef DEBUG
              mflux_out2(iblk) = mflux_out2(iblk) + mflux
#endif

            end if

          end do
        end do

      end if

C----------------------------------------------------------------------
C SCALE VELOCITIES: SCALE THE INTERMEDIATE VELOCITY IF task=1;
C                   SCALE THE VELOCITY AT THE NEW TIME STEP IF task=2
C----------------------------------------------------------------------

      if (iva .eq. SV_U) then
#ifdef DEBUG
        call usr_sum_double (mflux_in2, mflux_in, numbl_total,
     &  REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
        red = mflux_in
        call MPI_ALLREDUCE(red, mflux_in, 1, REAL_TYPE_MPI,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      endif

      if (iva .eq. SV_U) then
#ifdef DEBUG
        call usr_sum_double (mflux_out2, mflux_out, numbl_total,
     &  REAL_TYPE_MPI, MPI_COMM_WORLD, 0, bl2pe, pe, .true., ierr)
#else
        red = mflux_out
        call MPI_ALLREDUCE(red, mflux_out, 1, REAL_TYPE_MPI,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      endif

C      if (iva .eq. 3) then
      if (iva .eq. SV_W) then
C Unlike calculation of mflux_in and mflux_out, MPI operation of 
C mchg_rate is within subroutine mass_change instead of here
        if (task .eq. 1) then
          call mass_change(rho,rhonp1,mchg_rate,rdidjp,dkw)
        endif

        do iva1 = SV_U,SV_W
          do iblk = 1,numbl
            do iside = 1,6

              if (bl_bc(iside,iblk) .eq. BCT_CONVECTIVE) then

                nip1 = mg_nijkp1(1,iblk,1)
                njp1 = mg_nijkp1(2,iblk,1)
                nkp1 = mg_nijkp1(3,iblk,1)

               if (task .eq. 1) then
                  scale_factor = (mflux_in - mchg_rate)/mflux_out

               endif


C NOTE: POINTER p3 IS GIVEN IN A DIFFERENT WAY FROM THE ABOVE

                p3 = 1 + mg_off(iblk,1) + bl_tlen_sg_ijk*(iva1-1)

                if (task .eq. 1) then
                  call conv_exitbc_scale(iva1,iside,iblk,nip1,njp1,nkp1,
     &                 u_int(p3),scale_factor)
                end if

              end if

            end do
          end do
        end do

        do iblk = 1,numbl
          do iside = 1,6
            mfluxbl_out2(iside,iblk) = 0.
          end do
        end do

      end if
 
      
      return
      end





