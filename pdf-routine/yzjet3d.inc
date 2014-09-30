!        Mesh Size
         IMPLICIT NONE
         INTEGER :: ierr, nprocs, myrank, NQ
         INTEGER L1, L2, L3, Q1, Q2, Q3, M1, M2, M3, N1, N2, N3
         PARAMETER (L1=1, L2=1, L3=1) ! numbers of processors
         PARAMETER (Q1=512, Q2=128, Q3=128) ! mesh
         PARAMETER (M1=Q1/L1, M2=Q2/L2, M3=Q3/L3)
         PARAMETER (N1=M1+1, N2=M2+1, N3=M3+1)               
         INTEGER O1, O2, O3
         REAL Offset_x, Offset_y, Offset_z

         COMMON / GENTPVM / O1, O2, O3, Offset_x, Offset_y, Offset_z,myrank,ierr,nprocs,NQ 
!----------------------------------------------------------
!        U1 = u velocity
!        U2 = v velocity
!        U3 = w velocity
!        P1 = pressure
!        H1, H2, H3, R1 = historic terms
!        Wn  = work fields
         REAL U1(0:N1, 0:N2, 0:N3),                                    &
              U2(0:N1, 0:N2, 0:N3),                                    &
              U3(0:N1, 0:N2, 0:N3),                                    &
              H1(0:N1, 0:N2, 0:N3),                                    &
              H2(0:N1, 0:N2, 0:N3),                                    &
              H3(0:N1, 0:N2, 0:N3),                                    &
              R1(0:N1, 0:N2, 0:N3),                                    &
              T1(0:N1, 0:N2, 0:N3),                                    &
              pass1(0:N1, 0:N2, 0:N3),                                 &
              rNO(0:N1, 0:N2, 0:N3),                                    &
              rO3(0:N1, 0:N2, 0:N3),                                    &
              rNO2(0:N1, 0:N2, 0:N3),                                   &
              rO2(0:N1, 0:N2, 0:N3),                                    &
              Wa(0:N1, 0:N2, 0:N3),                                    &
              Wb(0:N1, 0:N2, 0:N3),                                    &
              Wc(0:N1, 0:N2, 0:N3),                                    &
              Wd(0:N1, 0:N2, 0:N3),            &
              We(0:N1, 0:N2, 0:N3),            &
              Wi(0:N1, 0:N2, 0:N3),            &
              Wg(0:N1, 0:N2, 0:N3),            &
              Wh(0:N1, 0:N2, 0:N3)

              

!        Field for large eddy simulation
!        E = eddy viscosity
         REAL E(0:N1, 0:N2, 0:N3)
         LOGICAL LES

         COMMON / FIELDS / U1, U2, U3, H1, H2, H3, Wa, Wb, Wc, Wd,     &
                           E, LES, rNO, rO3, rNO2, rO2,  &
                           R1, T1, pass1
!----------------------------------------------------------
!        SGS Models
         INTEGER None, Smagorinsky, F2, Structure, Germano
         PARAMETER (None=0, Smagorinsky=1, F2=2, Structure=3, Germano=4)
         INTEGER  Model, Model_Dim, Test_Fil
         REAL     Smag_Const_Sq
         LOGICAL  Low_Re_Mod
         REAL Rdy_T(1:M2), Rdy_B(1:M2), Rdx_L(1:M1), Rdx_R(1:M1)

         COMMON / MODELS / Model, Smag_Const_Sq, Low_Re_Mod, &
                           Model_Dim, Rdy_T, Rdy_B, Rdx_L, Rdx_R, &
                           Test_Fil
!----------------------------------------------------------        
!        Simulation control, time-stepping, etc
         LOGICAL Abort, New_Run
         INTEGER Last, Step_No, Debug, Count,statstart
         REAL    Time, Delta_t, dtover2, twoverdt, Alpha, Beta
         CHARACTER*30 Abort_Message
         REAL Unset
         PARAMETER (Unset = -1.11111111E11)

         COMMON / STEPS / Abort, New_Run, Last, Step_No, Debug, Time,  &
                          Delta_t, dtover2, twoverdt, Alpha, Beta, Count,&
                          statstart
         COMMON /ABORT /  Abort_Message
!----------------------------------------------------------
!        Diffusivities etc
         REAL Nu, NuonCv, Lambda, KonRhoCv, Inv_R       
         COMMON / DIFFUSE / Nu, NuonCv, Lambda, KonRhoCv, Inv_R
!----------------------------------------------------------
!        Global forces
         Real Force_x, Force_y, Force_z

         COMMON / FORCES / Force_x, Force_y, Force_z
!----------------------------------------------------------
!        Pressure solver & ffts
         REAL epsilon, accuracy
         REAL ksqu(Q3), Re(0:Q3/2), Im(0:Q3/2)  ! CPVM Qn
         INTEGER Max_Iterations
         INTEGER posf(Q3), negf(Q3)             ! CPVM Qn
         PARAMETER (Max_Iterations = 20)

         COMMON / PRESSR / epsilon, Accuracy, ksqu, posf, negf, Re, Im
!----------------------------------------------------------
!        Output control, statistics, etc
         INTEGER Max_Outs, Max_Planes,  Max_Traces
         PARAMETER (Max_Outs = 16, Max_Planes = 1, Max_Traces = 9)
         LOGICAL Out(1:Max_Outs)
         INTEGER XPl(1:Max_Planes),YPl(1:Max_Planes),ZPl(1:Max_Planes)
         INTEGER XTr(1:Max_Traces),YTr(1:Max_Traces),ZTr(1:Max_Traces)
         INTEGER X_Planes, Y_Planes, Z_Planes, Trace_Pts,              &
                 Dump_Freq, Stat_Freq, Moni_Freq, Fiel_Freq, Sour_Freq,&
                 Trac_Freq, No_Sam, Spec_Freq

         COMMON / OUTPUT / Out, XPl, YPl, ZPl, XTr, YTr, ZTr,          &
                           X_Planes, Y_Planes, Z_Planes, Trace_Pts,    &
                           Dump_Freq, Stat_Freq, Moni_Freq, Fiel_Freq, &
                           Sour_Freq,Trac_Freq, No_Sam, Spec_Freq                 

         REAL Courant_no, Viscous_no, CPUtot_n, CPUtot_o, CPUtot_d,    &
              CPUpr_n, CPUpr_o, CPUpr_d, SECOND

         COMMON / STATS / Courant_no, Viscous_no, CPUtot_n, CPUtot_o,  &
                          CPUtot_d, CPUpr_n, CPUpr_o, CPUpr_d 


!Added in extra stats wam3 8/11/04
         REAL u1_acc(1:M1,1:M2),u2_acc(1:M1,1:M2),u3_acc(1:M1,1:M2),   &
              p1_acc(1:M1,1:M2),T1_acc(1:M1,1:M2),                     &
              uu_acc(1:M1,1:M2),vv_acc(1:M1,1:M2),ww_acc(1:M1,1:M2),   &
              uv_acc(1:M1,1:M2),uw_acc(1:M1,1:M2),vw_acc(1:M1,1:M2),   &
              up_acc(1:M1,1:M2),vp_acc(1:M1,1:M2),wp_acc(1:M1,1:M2),   &
              pp_acc(1:M1,1:M2),tt_acc(1:M1,1:M2),T1_flu(1:M1,1:M2),   &
              cx_acc(1:M1,1:M2),cy_acc(1:M1,1:M2),cz_acc(1:M1,1:M2),   &
              uu_flu(1:M1,1:M2),vv_flu(1:M1,1:M2),ww_flu(1:M1,1:M2),   &
              uv_flu(1:M1,1:M2),uw_flu(1:M1,1:M2),vw_flu(1:M1,1:M2),   &
              pass1_acc(1:M1,1:M2),pp1_acc(1:M1,1:M2),                  &
              pass1_flu(1:M2,1:M2),rNO_acc(1:M1,1:M2),                &
              ppNO_acc(1:M1,1:M2),rNO_flu(1:M1,1:M2),                  &
              rO3_acc(1:M1,1:M2),ppO3_acc(1:M1,1:M2),rO3_flu(1:M1,1:M2),&
              rNO2_acc(1:M1,1:M2),ppNO2_acc(1:M1,1:M2),rNO2_flu(1:M1,1:M2),&
              rO2_acc(1:M1,1:M2),ppO2_acc(1:M1,1:M2),rO2_flu(1:M1,1:M2)
                 

         COMMON / ACCUM / u1_acc,u2_acc,u3_acc,p1_acc, &
                          uu_acc,vv_acc,ww_acc,uv_acc,uw_acc,vw_acc,   &
                       up_acc,vp_acc,wp_acc,pp_acc,cx_acc,cy_acc,cz_acc, &
                       uu_flu,vv_flu,ww_flu,uv_flu,uw_flu,vw_flu,tt_acc, &
                       T1_flu, T1_acc, pass1_acc, pp1_acc,pass1_flu,     &
                       rNO_acc,ppNO_acc,rNO_flu,rO3_acc,ppO3_acc,rO3_flu,&
                       rNO2_acc,ppNO2_acc,rNO2_flu,rO2_acc,ppO2_acc,     &
                       rO2_flu
!----------------------------------------------------------
!        Boundary conditions
         INTEGER Max_Boundaries
         PARAMETER (Max_Boundaries = 7)
         INTEGER Surface, Normal, Fixed, Gradient, Inflow
         INTEGER Cyclic, Convective, PVM_talk, jsplit
         PARAMETER (Surface = 1, Normal = -1, Fixed = 0) 
         PARAMETER (Gradient = 2, Inflow = 3, Cyclic = 4 )
         PARAMETER (Convective = 5, PVM_talk = 100)
         INTEGER Period(3), Boundaries
         INTEGER  Ilo(1:Max_Boundaries),                               &
                  Ihi(1:Max_Boundaries),                               &
                  Jlo(1:Max_Boundaries),                               &
                  Jhi(1:Max_Boundaries),                               &
                  Klo(1:Max_Boundaries),                               &
                  Khi(1:Max_Boundaries)
         REAL     bc(1:Max_Boundaries, 1:5),                           &
                  gr(1:Max_Boundaries, 1:5)                            
         INTEGER  sl(1:Max_Boundaries, 1:5),                           &
                  Type(1:Max_Boundaries, 1:5)
         INTEGER  Orient(1:Max_Boundaries, 1:3)

!     ####### not used at present
!         INTEGER  XSlices, XSliced
!         PARAMETER (XSlices=1, XSliced=1)
!         REAL XSlice(1:XSlices, 1:XSliced, 0:N2, 0:N3)
!         INTEGER YSlices, YSliced
!         PARAMETER (YSlices=1, YSliced=1)
!         REAL YSlice(1:YSlices, 1:YSliced, 0:N1, 0:N3)
!         INTEGER ZSlices, ZSliced
!         PARAMETER (ZSlices=1, ZSliced=1)
!         REAL ZSlice(1:ZSlices, 1:ZSliced, 0:N1, 0:N2)
!     #######

         COMMON / BOUNDS / Period, Boundaries,jsplit,                  &
                           Ilo, Ihi, Jlo, Jhi, Klo, Khi,               &
                           bc, gr, sl, Type, Orient                   ! &
!                           ,XSlice, YSlice, ZSlice                    

!     planes for inflow boundary conditions & profile for convective 
!     outflowc
!
!     #changed
      REAL U1jet_in(0:N2),  U2jet_in(0:N2), T_in(0:N2),              &
      U2en_top(1:M1), U2en_bot(1:M1),                                &
      U1_ran(0:N2,0:N3), U2_ran(0:N2,0:N3), U3_ran(0:N2,0:N3)
!     #changed
!
      COMMON /INFLOWC/ U1jet_in, U2jet_in, U2en_top, U2en_bot,       &
      U1_ran, U2_ran, U3_ran, T_in
!
!     next 2 lines not used
      REAL SUMOUTBL,U1bl (0:N2)     !constants for a simulation

      COMMON /OUTFLOWC/ U1bl, SUMOUTBL
!----------------------------------------------------------
!        Geometry
         REAL Delta_x, Size_x, Origin_x, Rdx
         REAL  x(0:N1+1), xsg(-1:N1)
         REAL dx(0:N1) 
         REAL cx0(0:N1), cx1(0:N1),                                    &
              xc0(0:N1), xc1(0:N1),                                    &
              dx0(0:N1), dx1(0:N1),                                    &
              tx0(0:N1), tx1(0:N1),                                    &
              xv(0:N1), tx(0:N1), dyu(0:N2), dxv(0:N1),                 &
              dxr(0:N1) 

         COMMON / XMESH / Delta_x, Size_x, Origin_x, Rdx, x, dx,       &
                          cx0, cx1, xc0, xc1, dx0, dx1, xsg, tx0,     &
                          tx1, xv, tx

         REAL Delta_y, Size_y, Origin_y, Rdy
         REAL  y(0:N2+1), ysg(-1:N2)
         REAL dy(0:N2)
         REAL cy0(0:N2), cy1(0:N2),                                    &
              yc0(0:N2), yc1(0:N2),                                    &
              dy0(0:N2), dy1(0:N2),                                    &
              ty0(0:N2), ty1(0:N2),                                    &
              yu(0:N2), ty(0:N2) 

         COMMON / YMESH / Delta_y, Size_y, Origin_y, Rdy, y, dy,       &
                          cy0,cy1, yc0,yc1, dy0,dy1, ty0, ty1, ysg,    &
                          yu, ty, dxv, dyu, dxr
         
         REAL Delta_z, Size_z, Origin_z, Rdz
         REAL z(0:N3+1), zsg(-1:N3)
         REAL dz(0:N3), tz(0:N2)  

         COMMON / ZMESH / Delta_z, Size_z, Origin_z, Rdz, z, dz, tz,  &
                          zsg 
!----------------------------------------------------------
!     Zebra solver
      INTEGER tall, mixt, wide
      PARAMETER (tall=0, mixt=1, wide=2)
      REAL roundoff 
      PARAMETER (roundoff=2.0E-8)

      INTEGER NxA, NyA, NxB, NyB, NxC, NyC
      PARAMETER (NxA=M1   +1, NyA=M2   +1)
      PARAMETER (NxB=M1/2 +1, NyB=M2/2 +1)
      PARAMETER (NxC=M1/4 +1, NyC=M2/4 +1)
!     etc ...
   
!     P1 is the pressure field. Wa, Wb, Wc are used for R1, B1, C1
      REAL P1(0:NxA,0:NyA,0:N3), G1(0:NxA,0:NyA)
      REAL Est1(0:NxA,0:NyA), Wst1(0:NxA,0:NyA), Nor1(0:NxA,0:NyA),    &
           Sou1(0:NxA,0:NyA), A1(0:NxA,0:NyA)
      REAL P2(0:NxB,0:NyB,0:N3), R2(0:NxB,0:NyB,0:N3),                 &
           B2(0:NxB,0:NyB,0:N3), C2(0:NxB,0:NyB,0:N3), G2(0:NxB,0:NyB)
      REAL Est2(0:NxB,0:NyB), Wst2(0:NxB,0:NyB), Nor2(0:NxB,0:NyB),    &
           Sou2(0:NxB,0:NyB), A2(0:NxB,0:NyB)
      REAL P3(0:NxC,0:NyC,0:N3), R3(0:NxC,0:NyC,0:N3),                 &
           B3(0:NxC,0:NyC,0:N3), C3(0:NxC,0:NyC,0:N3), G3(0:NxC,0:NyC)
      REAL Est3(0:NxC,0:NyC), Wst3(0:NxC,0:NyC), Nor3(0:NxC,0:NyC),    &
           Sou3(0:NxC,0:NyC), A3(0:NxC,0:NyC)

      INTEGER cells
      REAL max_residual, threshold, ARmax, ARmin

      COMMON /zebra/ max_residual, threshold, cells, ARmax, ARmin,     &
             P1, P2, P3, R2, R3, B2, B3, C2, C3, G1, G2, G3,           &
             Est1, Est2, Est3, Wst1, Wst2, Wst3, Nor1, Nor2, Nor3,     &
             Sou1, Sou2, Sou3, A1, A2, A3
!---------------------------------------------------------------------§
!        Eulerian PDF arrays
		 REAL :: pasf1(0:N1, 0:N2, 0:N3), nfields, NOin, O3in
         real :: rrate1, rrate2, rrate3, rrate4
         integer :: NOstream, passstream
         real :: NOf(0:N1, 0:N2, 0:N3), O3f(0:N1, 0:N2, 0:N3),    &
                 NO2f(0:N1, 0:N2, 0:N3), O2f(0:N1, 0:N2, 0:N3)
	    REAL :: areawf(0:N1,0:N2,0:N3),areaef(0:N1,0:N2,0:N3),   &
		        areasf(0:N1,0:N2,0:N3),areanf(0:N1,0:N2,0:N3),    &
				arealf(0:N1,0:N2,0:N3),arearf(0:N1,0:N2,0:N3)

		 REAL :: pass1M, sclam, scturb

         COMMON / EULERPDF / pasf1,pass1M,sclam,scturb,     &
				areawf,areaef,areanf,areasf,arealf,arearf,nfields, &
                NOin, O3in, NOstream, NOf, O3f, NO2f, O2f,         &
                rrate1,rrate2,rrate3,rrate4, passstream
!------------------------- END ----------------------------
