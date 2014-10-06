!----------------------------------------------------------------------
! gentles: free shear layer simulation code.
!----------------------------------------------------------------------
! Addition by Andrew McMullan
!    1.  Passive scalar routine (SUBROUTINE TEMPERATURE).
!    2.  Eulerian Stochastic Fields method
!----------------------------------------------------------------------
!
!
      PROGRAM Gentles
!     ---------------

!      (c) Peter Voke 1992, 1994
!      (c) Andrew McMullan, 2004-2007
!      
       INCLUDE 'jet3d.inc'
       INCLUDE 'mpif.h'
       INTEGER, DIMENSION(:), ALLOCATABLE :: seed
       integer :: numb
       double precision :: looptime, tlooptime
       double precision :: t_loop_start, t_loop_end

       CALL MPI_INIT(ierr)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
       WRITE(*,*) 'my rank is', myrank
       write(*,*) 'there are', nprocs, 'cpus'


      CALL RANDOM_SEED(SIZE=numb)
      ALLOCATE(seed(numb))
! change this.
      OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
      READ(89) seed
      CLOSE(89)
      CALL RANDOM_SEED(put=seed)

       if(myrank .eq. 7) write(*,*) 'hello!'
!       call random_seed(myrank)

       CALL Start

No_Sam=1
if(myrank .eq. 0) then
OPEN(UNIT=99,FILE='singplanviz.dat')
OPEN(UNIT=90,FILE='singplanscalviz.dat')
OPEN(UNIT=65,FILE='spanavaevort.dat')
OPEN(UNIT=91,FILE='spanaveviz.dat')
OPEN(UNIT=92,FILE='yz1.dat')
OPEN(UNIT=93,FILE='yz2.dat')
OPEN(UNIT=94,FILE='yz3.dat')
OPEN(UNIT=95,FILE='yz4.dat')
OPEN(UNIT=96,FILE='yz5.dat')
OPEN(UNIT=97,FILE='yz6.dat')
endif


        do while(abort .eqv. .false.)
        t_loop_start = MPI_WTIME()

         CALL Step

        t_loop_end = MPI_WTIME()
       looptime = t_loop_end-t_loop_start
       call MPI_REDUCE(looptime, tlooptime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)


       if(myrank .eq. 0)then
       write(*,'(''*Info: Time used: '',g15.8, '' sec. for time step '',i10)') tlooptime,Step_No-1
       end if
        enddo

       CALL Finish
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_FINALIZE(ierr)

       STOP       

      END

!----------------------------------------------------------------------

       SUBROUTINE Start
!      ---------------- 
         INCLUDE 'jet3d.inc'
!PVM

         CALL Assign_Files
         CALL Set_PVM
         CALL Defaults
         CALL Data_Read
         CALL fft_Table
         call face_areas
         CALL MG_Table
         IF ( New_Run ) THEN 
           CALL Set_New
         ELSE
           CALL Restart
         ENDIF
         CALL Read_inflow
!# not used         CALL Read_outflow
!         CALL Monitor

         RETURN
       END

!----------------------------------------------------------------------
       
       SUBROUTINE Set_PVM
!      ------------------
!PVM   Sets up the PVM variables

         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         INTEGER ip

!        Numbers of processors in 3D array
         Lproc(1) = L1
         Lproc(2) = L2
         Lproc(3) = L3

!        Find coordinates of processor and set up array me(1:3)
!        where (me(1),me(2),me(3)) = (0,0,0) -> (L1-1,L2-1,L3-1)
!        for a 3D array of L1 * L2 * L3 processors.

!        Set values of O1, O2, O3 = (M1*me(1),M2*me(2),M3*me(3))

!        Find tids of all colinear processors (3D cyclic please)
!        and put them into Xtid(-16:16), Ytid(-16:16), Ztid(-16:16)
!        Xtid(-1) = south, Xtid(1) = north
!        Ytid(-1) = east,  Ytid(1) = west
!        Ztid(-1) = below, Ztid(1) = above.
!        Ztid(-k) is tid of kth processor down.
!        Ztid(+k) is tid of kth processor up.
!        Xtid(0) = Ytid(0) = Ztid(0) = mytid

!        defaults for one-processor system follow:
         O1 = 0
         O2 = 0
         O3 = 0
         me(1) = 0
         me(2) = 0
         me(3) = 0
         DO ip = -16, 16
           Xtid(ip) = 0
           Ytid(ip) = 0
           Ztid(ip) = 0
         ENDDO

         RETURN
       END

!----------------------------------------------------------------------
      SUBROUTINE Face_areas
!     ----------------------


! Set face areas for PDF calculation - save multiple computations
	  
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'	  
	     INTEGER :: i,j,k
	  
	  
	DO K=1, M3
		DO J=1, M2
			DO I=1, M1
			areawf(i,j,k) = dz(k) * dyu(j)
			areaef(i,j,k) = dz(k) * dyu(j)
			areasf(i,j,k) = dxv(i) * dz(k)
			areanf(i,j,k) = dxv(i) * dz(k)
			arealf(i,j,k) = dxv(i) * dyu(j)
			arearf(i,j,k) = dxv(i) * dyu(j)
			ENDDO
		ENDDO
	ENDDO
	
	
	RETURN
	END
!----------------------------------------------------------------------
       SUBROUTINE Defaults
!      -------------------
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         INTEGER bo, pl, tr, which
 
         Nu = 0.001
         Last = 10
         Delta_t = 0.01
         New_Run = .TRUE.
         Count = 1
         Debug = 0
         Model = 0
         Low_Re_Mod = .FALSE.
         epsilon = 1.0E-4
         accuracy = 1.0E-6 !allowed maximum of local divergence (incomp.)
         Origin_x = 0.0
         Origin_y = 0.0
         Origin_z = 0.0
         Force_x = 0.0
         Force_y = 0.0
         Force_z = 0.0
         Stat_Freq = 10**4
         Moni_Freq = 10**4
         Fiel_Freq = 10**4
         Dump_Freq = 10**4      
         Sour_Freq = 10**4
         Trac_Freq = 10**4

! Increment bases
         Boundaries = 0   
         X_Planes = 0
         Y_Planes = 0
         Z_Planes = 0
         Trace_Pts= 0
 
! Broadcast array assignments
        DO bo = 1, max_boundaries
         Out(bo) = .FALSE.
         Ilo(bo) = 1
         Ihi(bo) = M1
         Jlo(bo) = 1
         Jhi(bo) = M2
         Klo(bo) = 1
         Khi(bo) = M3
         DO which = 1, 5
!PVM      default is PVM-Talk (even for pressure)
          bc(bo,which) = 0.0       
          gr(bo,which) = Unset
          Type(bo,which) = PVM_talk
          sl(bo,which) = 0
         ENDDO 
         Orient(bo,1) = 0
         Orient(bo,2) = 0
         Orient(bo,3) = 0
        ENDDO ! bo
!PVM Pressure Neumann conditions on t-r-u-e external boundaries
!    (Needed for internal boundaries too)
        IF (me(1) .EQ. 0) THEN
          Type(1,5) = Gradient
          gr(1,5) = 0.0  
        ENDIF
        IF (me(1) .EQ. L1-1) THEN
          Type(2,5) = Gradient
          gr(2,5) = 0.0  
        ENDIF
        IF (me(2) .EQ. 0) THEN
          Type(3,5) = Gradient
          gr(3,5) = 0.0  
        ENDIF
        IF (me(2) .EQ. L2-1) THEN
          Type(4,5) = Gradient
          gr(4,5) = 0.0  
        ENDIF
        IF (me(3) .EQ. 0) THEN
          Type(5,5) = Gradient
          gr(5,5) = 0.0  
        ENDIF
        IF (me(3) .EQ. L3-1) THEN
          Type(6,5) = Gradient
          gr(6,5) = 0.0  
        ENDIF

!PVM     default external boundaries
         Ilo(1) = 1   ! bo=1 is left
         Ihi(2) = M1  ! bo=2 is right
         Jlo(3) = 1   ! bo=3 is bottom
         Jhi(4) = M2  ! bo=4 is top
         Klo(5) = 1   ! bo=5 is near
         Khi(6) = M3  ! bo=6 is far
         Orient(1,1) = +1
         Orient(2,1) = -1
         Orient(3,2) = +1
         Orient(4,2) = -1
         Orient(5,3) = +1
         Orient(6,3) = -1

        DO pl = 1, Max_Planes
         XPl(pl) = 0
         YPl(pl) = 0
         ZPl(pl) = 0
        ENDDO
!        DO tr = 1, Max_Traces
!         XTr(tr) = 0
!         YTr(tr) = 0
!         ZTr(tr) = 0
!        ENDDO
!PVM     Qn-1
         Period(1) = Q1-1
         Period(2) = Q2-1
         Period(3) = Q3-1
!         CPUtot_o=SECOND()   !get cpu time
         CPUpr_d=0.0
 
         RETURN
 
       END
 
!----------------------------------------------------------------------

       SUBROUTINE Data_Read
!      -------------------- 
                              
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         CHARACTER*30 String30, Temp30
         CHARACTER* 6 String, Temp
         REAL Value
         INTEGER Round, bo, face, i, j, k, which
         LOGICAL Intersects

         WRITE (*,*) 'DATA_READ subroutine...' 
         Intersects = .TRUE.
         face = 1 ! prevents an error
         OPEN ( 10, FILE='jet3d.dat', STATUS='OLD' )

 10      CONTINUE
           READ (10, '(F20.10, A)', END=20 ) Value, String30
           IF ( ABS(Value) .GT. 0.0 ) THEN
           if(myrank .eq. 0) then
             PRINT *, String30, '  =  ', Value
           endif 
           ELSE
           if (myrank .eq. 0) then
             PRINT *, String30
           endif
           ENDIF
           Round = INT( Value+0.1 )
           
           IF (String30 .EQ. ' ') String30 = '*'
           DO WHILE ( String30(1:1) .EQ. ' ' )
             Temp30 = String30(2:30) // ' '
             String30 = Temp30
           ENDDO
           i = INDEX( String30,' ' )
           IF (i .GT. 0) THEN 
             Temp30 = String30(1:i-1)
             String30 = Temp30
           ENDIF
           String = String30(1:6)
           DO i = 1, 6
             k = ICHAR(String(i:i))
             IF (k .GT. 95) THEN
               k = k-32
               Temp = String(1:i-1) // CHAR(k) // String(i+1:6)
               String = Temp
             ENDIF
           ENDDO

             IF ( String .EQ. 'VISCOS') THEN
               Nu = Value 
               Lambda = -Nu * 2. / 3.    
             ENDIF
             IF ( String .EQ. 'LAST' )  Last = Round
             IF ( String .EQ. 'TIME' )  Delta_t = Value
             IF ( String .EQ. 'NEW' )   New_Run = .TRUE.
             IF ( String .EQ. 'OLD' )   New_Run = .FALSE.
             IF ( String .EQ. 'DUMP' )  Dump_Freq = Round
             IF ( String .EQ. 'DEBUG' ) Debug = Round
             IF ( String .EQ. 'DIRECT') THEN
               LES = .FALSE.
               Model = None                 
             ENDIF
             IF ( String .EQ. 'SMAGOR') THEN
               LES = .TRUE.
               Model = Smagorinsky
               Smag_Const_Sq = Value**2     
             ENDIF
             IF ( String .EQ. 'GERM') THEN
               LES = .TRUE.
               Model = Germano 
               Test_Fil = INT( Value )
             ENDIF
             IF ( String .EQ. 'STRUCT' ) THEN
               LES = .TRUE.
               Model = Structure
               Model_Dim = INT( Value )
             ENDIF
             IF ( String .EQ. 'F2' ) THEN
               LES = .TRUE.
               Model = F2
             ENDIF
             IF ( String .EQ. 'LOW' )   Low_Re_Mod = .TRUE.
             IF ( STRING .EQ. 'pass1M') pass1M = Value
             IF ( STRING .EQ. 'NOIN') NOin = Value
             IF ( STRING .EQ. 'O3IN') O3in = Value
             IF ( STRING .EQ. 'NOSTRE') NOstream = Round
             IF ( STRING .EQ. 'PASSST') passstream = Round
             IF ( String .EQ. 'RRATE1' )  rrate1 = Value
             IF ( String .EQ. 'RRATE2' )  rrate2 = Value
             IF ( String .EQ. 'RRATE3' )  rrate3 = Value
             IF ( String .EQ. 'RRATE4' )  rrate4 = Value
             IF ( String .EQ. 'SCTURB' )  Scturb = Value
             IF ( String .EQ. 'SCLAM' )  Sclam = Value
             IF ( String .EQ. 'ACCURA') epsilon = Value
             IF ( String .EQ. 'SPLITT' )  jsplit = Round
             IF ( String .EQ. 'NFIELD' )  nfields = Round
             IF ( String .EQ. 'STATST' )  statstart = Round
             IF ( String .EQ. 'DX' )    CALL X_Table ( 0, Value )
             IF ( String .EQ. 'DY' )    CALL Y_Table ( 0, Value )  
             IF ( String .EQ. 'DZ' )    CALL Z_Table ( 0, Value )
!PVM next 3 lines / Qn
             IF ( String .EQ. 'X-SIZE') CALL X_Table ( 0, Value / Q1 )
             IF ( String .EQ. 'Y-SIZE') CALL Y_Table ( 0, Value / Q2 )
             IF ( String .EQ. 'Z-SIZE') CALL Z_Table ( 0, Value / Q3 )
             IF ( String .EQ. 'X-MESH') CALL X_Table ( Round, 0.0 )
             IF ( String .EQ. 'Y-MESH') CALL Y_Table ( Round, 0.0 )
             IF ( String .EQ. 'Z-MESH') CALL Z_Table ( Round, 0.0 )
             IF ( String .EQ. 'X-SHIF') THEN
               Origin_x = value
               CALL X_Table ( 0, Delta_x )  
             ENDIF
             IF ( String .EQ. 'Y-SHIF') THEN
               Origin_y = value
               CALL Y_Table ( 0, Delta_y )
             ENDIF
             IF ( String .EQ. 'Z-SHIF') THEN
               Origin_z = value
               CALL Y_Table ( 0, Delta_y )
             ENDIF
             IF ( String .EQ. 'X-FORC') Force_x = Value
             IF ( String .EQ. 'Y-FORC') Force_y = Value
             IF ( String .EQ. 'Z-FORC') Force_z = Value
! boundary conditions
             IF ( String .EQ. 'BOUNDA') THEN
               bo = Boundaries + 1
               IF (bo .NE. Round .AND. Round .NE. 0) THEN
                 PRINT *, ' Warning, boundary ', Round,' out of order.'
               ENDIF
               IF (bo .GT. 6 .AND. .NOT. Intersects) bo = bo - 1
!PVM               as last boundary did not intersect
               Boundaries = bo              
               Intersects = .TRUE. ! assume next boundary intersects
             ENDIF

!PVM Lots of PVM adjustments follow.
!PVM Logical variable Intersects is used to flag intersection of 
!PVM a boundary with the processor subdomain.
!PVM Boundaries 1 to 6 are interprocessor talk boundaries by default,
!PVM but may be overwritten. 
             IF ( String .EQ. 'X+' ) THEN
               i = Round - O1
               IF (i .GE. 1 .AND. i .LE. M1 ) THEN
                 Ilo(bo) = i
                 Ihi(bo) = i
                 Orient(bo,1) = +1 
                 face = 1
                 Type(bo,5) = Gradient ! sets p.b.c. 
                 gr(bo,5) = 0.0  
               ELSE
                 Intersects = .FALSE. ! boundary doesn't intersect
               ENDIF
             ENDIF
             IF ( String .EQ. 'X-' ) THEN
               i = Round - O1
               IF (i .GE. 1 .AND. i .LE. M1) THEN
                 Ilo(bo) = i
                 Ihi(bo) = i
                 Orient(bo,1) = -1 
                 face = 1
                 Type(bo,5) = Gradient ! sets p.b.c. 
                 gr(bo,5) = 0.0  
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Y+' ) THEN
               j = Round - O2
               IF (j .GE. 1 .AND. j .LE. M2) THEN
                 Jlo(bo) = Round - O2
                 Jhi(bo) = Round - O2
                 Orient(bo,2) = +1 
                 face = 2
                 Type(bo,5) = Gradient ! sets p.b.c. 
                 gr(bo,5) = 0.0  
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Y-' ) THEN
               j = Round - O2
               IF (j .GE. 1 .AND. j .LE. M2) THEN  !CHECK THIS
                 Jlo(bo) = j
                 Jhi(bo) = j
                 Orient(bo,2) = -1 
                 face = 2
                 Type(bo,5) = Gradient ! sets p.b.c. 
                 gr(bo,5) = 0.0  
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Z+' ) THEN
               k = Round - O3
               IF (k .GE. 1 .AND. k .LE. M3) THEN
                 Klo(bo) = Round - O3
                 Khi(bo) = Round - O3
                 Orient(bo,3) = +1 
                 face = 3
                 Type(bo,5) = Gradient ! sets p.b.c. 
                 gr(bo,5) = 0.0  
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Z-' ) THEN
               k = Round - O3
               IF (k .GE. 1 .AND. k .LE. M3) THEN
                 Klo(bo) = Round - O3
                 Khi(bo) = Round - O3
                 Orient(bo,3) = -1 
                 face = 3
                 Type(bo,5) = Gradient ! sets p.b.c. 
                 gr(bo,5) = 0.0  
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'X-LOW' .AND. Intersects ) THEN
               i = Round - O1
               IF (i .LE. M1) THEN
                 Ilo(bo) = i
               ELSE
                 Intersects = .FALSE. ! boundary doesn't intersect
               ENDIF
             ENDIF
             IF ( String .EQ. 'X-HI' .AND. Intersects ) THEN   
               i = Round - O1
               IF (i .GE. 0) THEN
                 Ihi(bo) = i
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Y-LOW' .AND. Intersects ) THEN  
               j = Round - O2
               IF (j .LE. M2) THEN
                 Jlo(bo) = j
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Y-HI' .AND. Intersects ) THEN  
               j = Round - O2
               IF (j .GE. 0) THEN
                 Jhi(bo) = j
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Z-LOW' .AND. Intersects ) THEN  
               k = Round - O3
               IF (k .LE. M3) THEN
                 Klo(bo) = k
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
             IF ( String .EQ. 'Z-HI' .AND. Intersects ) THEN
               k = Round - O3
               IF (i .GE. 0) THEN
                 Khi(bo) = k
               ELSE
                 Intersects = .FALSE.
               ENDIF
             ENDIF
! no slip, surface value
             IF ( String .EQ. 'U-SURF'.AND. Intersects ) THEN
               bc(bo,1) = Value
               Type(bo,1) = Surface
               IF (face .EQ. 1 ) Type(bo,1) = Normal
             ENDIF
             IF ( String .EQ. 'V-SURF'.AND. Intersects ) THEN
               bc(bo,2) = Value
               Type(bo,2) = Surface
               IF (face .EQ. 2 ) Type(bo,2) = Normal
             ENDIF
             IF ( String .EQ. 'W-SURF'.AND. Intersects ) THEN
               bc(bo,3) = Value
               Type(bo,3) = Surface
               IF (face .EQ. 3 ) Type(bo,3) = Normal
             ENDIF
             IF ( String .EQ. 'T-SURF'.AND. Intersects ) THEN
               bc(bo,4) = Value
               Type(bo,4) = Surface
             ENDIF
! slip, normal gradient * delta-n value
             IF ( String .EQ. 'U-GRAD'.AND. Intersects ) THEN
               gr(bo,1) = Value
               Type(bo,1) = Gradient
             ENDIF
             IF ( String .EQ. 'V-GRAD'.AND. Intersects ) THEN
               gr(bo,2) = Value
               Type(bo,2) = Gradient
             ENDIF
             IF ( String .EQ. 'W-GRAD'.AND. Intersects ) THEN
               gr(bo,3) = Value
               Type(bo,3) = Gradient
             ENDIF
             IF ( String .EQ. 'T-GRAD'.AND. Intersects ) THEN
               gr(bo,4) = Value
               Type(bo,4) = Gradient
             ENDIF
! sliced
             IF ( String .EQ. 'U-FILE'.AND. Intersects ) THEN
               sl(bo,1) = Round
               Type(bo,1) = Inflow
             ENDIF
             IF ( String .EQ. 'V-FILE'.AND. Intersects ) THEN
               sl(bo,2) = Round
               Type(bo,2) = Inflow
             ENDIF
             IF ( String .EQ. 'W-FILE'.AND. Intersects ) THEN
               sl(bo,3) = Round
               Type(bo,3) = Inflow
             ENDIF
             IF ( String .EQ. 'T-FILE'.AND. Intersects ) THEN
               sl(bo,4) = Round
               Type(bo,4) = Inflow
             ENDIF
! cyclic
             IF ( String .EQ. 'CYCLIC'.AND. Intersects                 &
                  .AND. Lproc(face) .EQ. 1) THEN
!PVM           If more than one processor for this dimension,
!              leave as a PVM_talk boundary (talk is cyclic)   
               Type(bo,1) = Cyclic
               Type(bo,2) = Cyclic
               Type(bo,3) = Cyclic
               Type(bo,5) = Cyclic
!              (as pressure is periodic too)
             ENDIF
             IF ( String .EQ. 'T-CYCL'.AND. Intersects                 &
                  .AND. Lproc(face) .EQ. 1) THEN
               Type(bo,4) = Cyclic
             ENDIF
! convective
             IF ( String .EQ. 'CONVEC'.AND. Intersects ) THEN
               Type(bo,1) = Convective
               Type(bo,2) = Convective
               Type(bo,3) = Convective
               Type(bo,4) = Convective
             ENDIF

! Statistics
             IF ( String .EQ. 'STATIS') Stat_Freq = Round
             IF ( String .EQ. 'MONITO') Moni_Freq = Round
             IF ( String .EQ. 'FIELDS') Fiel_Freq = Round
             IF ( String .EQ. 'ACOUST') Sour_Freq = Round
             IF ( String .EQ. 'SPECTR') Spec_Freq = Round
             IF ( String .EQ. 'TRACEP')  Trac_Freq = Round
             IF ( String .EQ. 'U' )      Out(1) = .TRUE.
             IF ( String .EQ. 'V' )      Out(2) = .TRUE.
             IF ( String .EQ. 'W' )      Out(3) = .TRUE.
             IF ( String .EQ. 'P' )      Out(4) = .TRUE.

             IF ( String .EQ. 'UU' )     Out(5) = .TRUE.
             IF ( String .EQ. 'VV' )     Out(6) = .TRUE.
             IF ( String .EQ. 'WW' )     Out(7) = .TRUE.
             IF ( String .EQ. 'PP' )     Out(8) = .TRUE.
                      
             IF ( String .EQ. 'UV' )     Out(9) = .TRUE.
             IF ( String .EQ. 'UW' )     Out(10) = .TRUE.
             IF ( String .EQ. 'VW' )     Out(11) = .TRUE.
             IF ( String .EQ. 'UP' )     Out(12) = .TRUE.
             IF ( String .EQ. 'VP' )     Out(13) = .TRUE.
             IF ( String .EQ. 'WP' )     Out(14) = .TRUE.
             IF ( String .EQ. 'CORI' )   Out(15) = .TRUE.   !Coriolis Acc.
             IF ( String .EQ. 'HELI' )   Out(16) = .TRUE.   !Helicity

             IF ( String .EQ. 'X-PLAN') THEN 
               IF (X_Planes .LT. Max_Planes ) THEN
                 X_Planes = X_Planes + 1
                 XPl(X_Planes) = Round
               ELSE
                 Abort_Message = ' Too many X planes.'
                 Abort = .TRUE. 
                 GOTO 15
               ENDIF
             ENDIF
             IF ( String .EQ. 'Y-PLAN') THEN
               IF (Y_Planes .LT. Max_Planes ) THEN
                 Y_Planes = Y_Planes + 1
                 YPl(Y_Planes) = Round
               ELSE
                 Abort_Message = ' Too many Y planes.'
                 Abort = .TRUE. 
                 GOTO 15
               ENDIF
             ENDIF
             IF ( String .EQ. 'Z-PLAN') THEN
               Z_Planes = Z_Planes + 1
               IF (Z_Planes .LT. Max_Planes ) THEN
                 Z_Planes = Z_Planes + 1
                 ZPl(Z_Planes) = Round
               ELSE
                 Abort_Message = 'Too many Z planes.'
                 Abort = .TRUE. 
                 GOTO 15
               ENDIF
             ENDIF
             IF ( String .EQ. 'TRACEP') THEN 
               IF (Trace_Pts .LT. Max_Traces) THEN
                 Trace_Pts = Trace_Pts + 1
               ELSE
                 Abort_Message = 'Too many trace points.'
                 Abort = .TRUE. 
                 GOTO 15
               ENDIF
             ENDIF
             IF ( String .EQ. 'X-TRAC') XTr(Trace_Pts) = Round
             IF ( String .EQ. 'Y-TRAC') YTr(Trace_Pts) = Round
             IF ( String .EQ. 'Z-TRAC') ZTr(Trace_Pts) = Round
             IF ( String .EQ. '*' ) String = '*' 
             IF ( String .EQ. ' ' ) String = ' '
!            do nothing: comment only
         IF ( String .NE. 'END' ) GOTO 10

 15      CONTINUE
         twoverdt = 1.0 / Delta_t
         dtover2  = Delta_t / 1.0
         Alpha  =  1.5 * Delta_t
         Beta   = -0.5 * Delta_t

!
         WRITE (*,*) myrank, '... DATA_READ completed!'

         if(myrank .eq. 0) then
         WRITE(6,*) ' Check Type(bo,which) '
         DO bo = 1,7
         DO which = 1,5
         WRITE(6,*) 'bo=',bo,'which=',which, 'Type(bo,which)=', &
         Type(bo,which)
         ENDDO
         ENDDO
!
         WRITE(6,*) ' Ilo, Ihi, Jlo, Jhi, Klo, Khi '
         DO bo = 1,7
         WRITE(6,*) 'bo=',bo,'Ilo=',Ilo(bo),'Ihi=',Ihi(bo)
         WRITE(6,*) 'bo=',bo,'Jlo=',Jlo(bo),'Jhi=',Jhi(bo)
         WRITE(6,*) 'bo=',bo,'Klo=',Klo(bo),'Khi=',Khi(bo)
         ENDDO
         endif !myrank
         RETURN

 20      PRINT *, '** Warning: abnormal termination of input file.'
         PRINT *, '** END missing from input file.'

         CLOSE (10)



         RETURN
 
       END

!----------------------------------------------------------------------

       SUBROUTINE Set_New
!      ------------------
         INCLUDE 'jet3d.inc'

         CALL Ini_U_P
         CALL Ini_W
         CALL Ini_H
         CALL Ini_S
         Step_No = 1
         Time = 0.0
         No_Sam = 1!    for statistics
         RETURN

       END

!----------------------------------------------------------------------

       SUBROUTINE Restart
!      ------------------
         INCLUDE 'jet3d.inc'
         INCLUDE 'mpif.h'
         INTEGER QN1, QN2, QN3, i, j, k
         LOGICAL LESQ
         REAL :: AO1(1:M1*M2*M3),AO2(1:M1*M2*M3),AO3(1:M1*M2*M3),AO4(1:M1*M2*M3)
         REAL :: AO5(1:M1*M2*M3),AO6(1:M1*M2*M3),AO7(1:M1*M2*M3)
         character*3 int_as_String
         character*12 all

         write(int_as_string, fmt='(I3.3)') myrank
         all='sfield'//int_as_string//'.3d'


!        1. read old velocity, pressure & history fields
         IF(myrank .EQ. 0) THEN
         OPEN ( UNIT=7, FILE='dump.3d', STATUS='OLD',               &
         FORM='UNFORMATTED' )
         REWIND (7)

         READ (7) Step_No, Time

         READ (7) QN1, QN2, QN3
         IF (QN1.NE.N1 .OR. QN2.NE.N2 .OR. QN3.NE.N3) THEN 
           Abort_message = 'Bad mesh size on restart (dump).'
           Abort = .TRUE. 
           RETURN
         ENDIF

         READ (7) LESQ
         IF ( LESQ .NEQV. LES ) THEN
           Abort_message = 'LES mismatched on restart.'
           Abort = .TRUE. 
           RETURN        
         ENDIF

         READ (7) (((U1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((U2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((U3(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((P1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((H1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((H2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((H3(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((T1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((pass1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((rNO(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((rO3(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((rNO2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         READ (7) (((rO2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         CLOSE (UNIT=7,STATUS='KEEP')
         ENDIF  !myrank

!Broadcast scalar to all cpus
IF(myrank .EQ. 0) THEN
DO I=0, N1
DO K=0, N3
DO J=0, N2
AO1( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = pass1(i,j,k)
ENDDO
ENDDO
ENDDO


ENDIF

!Broadcast it
CALL MPI_BCAST(AO1, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
!unpack it
IF(myrank .GT. 0) THEN
    DO I=0, N1
    DO K=0, N3
    DO J=0, N2
        pass1(i,j,k) = AO1( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
    ENDDO
    ENDDO
    ENDDO
ENDIF

!        2.  read stochastic fields
         open(unit=1000+myrank,form='unformatted')
         rewind(1000+myrank)
         read(1000+myrank) (((pasf1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         read(1000+myrank) (((NOf(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         read(1000+myrank) (((O3f(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         read(1000+myrank) (((NO2f(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         read(1000+myrank) (((O2f(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         close(unit=1000+myrank, status='keep')


!        2. read old statistics fields for accumulation

!      IF(myrank .EQ. 0) THEN
!OPEN ( UNIT=7, FILE='stat.3d')
!REWIND(7)
!READ(7,*) No_Sam
!DO I=1,M1
!DO J=1,M2
!READ(7,*) u1_acc(i,j),u2_acc(i,j),u3_acc(i,j),p1_acc(i,j),T1_acc(i,j)
!ENDDO
!ENDDO

!DO I=1,M1
!DO J=1,M2
!READ(7,*) uu_acc(i,j),vv_acc(i,j),ww_acc(i,j),pp_acc(i,j),vp_acc(i,j)
!ENDDO
!ENDDO
!DO I=1,M1
!DO J=1,M2
!READ(7,*) uv_acc(i,j),uw_acc(i,j),vw_acc(i,j),up_acc(i,j),wp_acc(i,j)
!ENDDO
!ENDDO
!DO I=1,M1
!DO J=1,M2
!READ(7,*) cx_acc(i,j),cy_acc(i,j),cz_acc(i,j),tt_acc(i,j),T1_flu(i,j)
!ENDDO
!ENDDO

!DO I=1,M1
!DO J=1,M2
!READ(7,*) uu_flu(i,j),vv_flu(i,j),ww_flu(i,j),uv_flu(i,j),uw_flu(i,j),vw_flu(i,j)
!ENDDO
!ENDDO


!DO I=1,M1
!DO J=1,M2
!READ(7,*) pass1_acc(i,j), pp1_acc(i,j),pass1_flu(i,j)
!ENDDO
!ENDDO

!DO I=1, M1
!DO J=1, M2
!READ(7,*) rNO_acc(i,j), ppNO_acc(i,j),rNO_flu(i,j)
!ENDDO
!ENDDO

!DO I=1, M1
!DO J=1, M2
!READ(7,*) rO3_acc(i,j), ppO3_acc(i,j),rO3_flu(i,j)
!ENDDO
!ENDDO

!DO I=1, M1
!DO J=1, M2
!READ(7,*) rNO2_acc(i,j), ppNO2_acc(i,j),rNO2_flu(i,j)
!ENDDO
!ENDDO

!DO I=1, M1
!DO J=1, M2
!READ(7,*) rO2_acc(i,j), ppO2_acc(i,j),rO2_flu(i,j)
!ENDDO
!ENDDO


!CLOSE (UNIT=7,STATUS='KEEP')
!endif !myrank

         RETURN

       END

!----------------------------------------------------------------------
 
       SUBROUTINE Step
!      ---------------
 
         INCLUDE 'jet3d.inc'
         INCLUDE 'mpif.h'
         LOGICAL acc_div
         REAL divergence

!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         IF(myrank .EQ. 0) THEN
         IF ( Debug .GT. 0 ) THEN
         WRITE (*,*) '############ Time Step= ', Step_No,' ############'
         ENDIF
        CALL Ini_Ran   !put random noises   
        CALL Velocity_bc (.TRUE.)  !set outflow bc as well      
         CALL SubGrid
         CALL Momentum
         CALL Velocity_bc (.FALSE.) !reset inflow & cyclic bc
        CALL Pressure
         CALL Temperature
         CALL Outflow_temp
         CALL Velocity_bc (.FALSE.) !reset inflow & cyclic bc
!         WRITE(*,*) 'Entering PDF calculation'
          ENDIF !myrank
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL Euler_PDF
         IF (MOD(Step_No, Dump_Freq) .EQ. 0 ) CALL Dump
	     IF(myrank .EQ. 0) THEN
         IF (MOD(Step_No, Moni_Freq) .EQ. 0 ) CALL Monitor
         IF (Step_No .ge. statstart ) then
         IF (MOD(Step_No, Stat_Freq) .EQ. 0 ) CALL Statistics
         IF (MOD(Step_No, Fiel_Freq) .EQ. 0 ) CALL Field_out
         IF (MOD(Step_No, Sour_Freq) .EQ. 0 ) CALL Acoustics
!         IF (MOD(Step_No, Trac_Freq) .EQ. 0 ) CALL Trace
         IF (MOD(Step_No, Spec_Freq) .EQ. 0 ) CALL Spectra
!        Continuity check
!         CALL Div_U (acc_div,divergence)
!         IF (acc_div .AND. Step_No .GT. 10 .AND. epsilon .GT. 1.0E-08)&
!         epsilon=epsilon/10.0
         endif !step_no
         endif !myrank
        
         Step_No = Step_No + 1
         Time = Time + Delta_t
         CALL Ini_W

         IF ( Step_No .GT. Last ) THEN
           Abort = .TRUE.
           Abort_Message = 'Normal termination'
         ENDIF !Step_No

         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


         RETURN
       END
!
!----------------------------------------------------------------------
!!!!!!
      SUBROUTINE Subgrid
!     ------------------
!     Smagorinsky's SGS Model and 
!     2nd Order Velocity Structure Function Model.
      INCLUDE 'jet3d.inc'
      INTEGER i, j, k
      REAL  Const, Emax, Ratio
      REAL Beta2, OneThird, TwoThird, OneSixth
      REAL dz_F, dz_B, Delta
      REAL Damp(1:N1, 1:N2, 1:N3), num(1:M1,1:M2), den(1:M1,1:M2)
      REAL sum, max,imax,jmax,kmax, Ev
      REAL :: frac, cc1,cc2,cc3, TIME1, TIME2

!
      Beta2 = 0.25
      OneThird = 1.0/3.0
      TwoThird = 2.0/3.0
      OneSixth = 1.0/6.0
!
      IF( Model .EQ. Smagorinsky ) THEN
!     SijSji = S11**2 + S22**2 + S33**2 + 
!            2( S12**2 + S23**2 + S13**2 )
!     SijSij is defined in the cell centre.

IF (Step_No .EQ. 1) THEN
DO i=2, M1
DO j=1, M2
DO k=1, M3
U1(i,j,k) = U1(1,j,k)
ENDDO
ENDDO
ENDDO
ENDIF    



      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     Wa = S11*S11 + S22*S22 + S33*S33 
      Wa(i,j,k) = ( (U1(i+1,j,k) - U1(i,j,k))*dx0(i) )**2 &
                + ( (U2(i,j+1,k) - U2(i,j,k))*dy0(j) )**2 &
                + ( (U3(i,j,k+1) - U3(i,j,k))*Rdz    )**2
      ENDDO
      ENDDO
      ENDDO


!

      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     S12 off center (at lower left corner node (i,j,k) )
      Wb(i,j,k) =  ( U1(i,j,k) - U1(i,j-1,k) )*dy1(j) &
                +  ( U2(i,j,k) - U2(i-1,j,k) )*dx1(i)   
      ENDDO
      ENDDO
      ENDDO

     
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     Interpolate S12 at the centre from 4 corner nodes 
!     Wa = Wa + 2*S12S12 
      Wa(i,j,k) = Wa(i,j,k)        &
      + 0.03125*( Wb(i,j,k)+Wb(i+1,j,k)+Wb(i,j+1,k)+Wb(i+1,j+1,k) )**2
      ENDDO
      ENDDO
      ENDDO


      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     S23 off center
      Wb(i,j,k) = ( U2(i,j,k) - U2(i,j,k-1) )*Rdz &
                + ( U3(i,j,k) - U3(i,j-1,k) )*dy1(j)  
      ENDDO
      ENDDO
      ENDDO


!

      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     Wa = Wa + 2*S23S23 at the cell center
      Wa(i,j,k) = Wa(i,j,k)          &
      + 0.03125*( Wb(i,j,k)+Wb(i,j+1,k)+Wb(i,j,k+1)+Wb(i,j+1,k+1) )**2
      ENDDO
      ENDDO
      ENDDO
    


      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     S13 off centre
      Wb(i,j,k) = ( U1(i,j,k) - U1(i,j,k-1) )*Rdz &
                + ( U3(i,j,k) - U3(i-1,j,k) )*dx1(i) 
      ENDDO
      ENDDO
      ENDDO



      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     Wa = Wa + 2*S13S13 at the cell center
      Wa(i,j,k) = Wa(i,j,k)         &
      + 0.03125*( Wb(i,j,k)+Wb(i+1,j,k)+Wb(i,j,k+1)+Wb(i+1,j,k+1) )**2
      ENDDO
      ENDDO
      ENDDO



!     Vs = ((Cs)**2) * (dx*dy*dz)**(2/3) * SQRT(2*SijSij)
!     E() = Nu + Vs

     
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      E(i,j,k) = Smag_Const_Sq * (dx(i)*dy(j)*dz(k) )**(TwoThird) &
              * SQRT( 2.0*Wa(i,j,k) ) + Nu 

      ENDDO
      ENDDO
      ENDDO


!     Change Nu to E()  in the momentum calculations!
!     Beta2=0.25 and TwoThird=2.0/3.0 

      Const = 1./(Beta2*Nu)
      IF(  Low_Re_Mod  ) THEN  ! Low Reynold number modification
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      E(i,j,k) = E(i,j,k)                   &
               - Beta2*Nu*(  1 - EXP( -( E(i,j,k)-Nu )*Const )  )
      ENDDO
      ENDDO
      ENDDO
      ENDIF
!
      ENDIF   ! endif of Smagorinsky's model

!  
!     Filtered 2nd Order Velocity Function Model:
!     ------------------------------------------
!     1-D model in the z-homogeneous direction;
!     3-D used only if accumulated time averages 
!     of velocity field are available.
!
    
      IF( Model .EQ. Structure )THEN
!

!     Interpolate velocity at the cell centre.
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      Wa(i,j,k) = ( U1(i,j,k) + U1(i+1,j,k) )*0.5 
      Wb(i,j,k) = ( U2(i,j,k) + U2(i,j+1,k) )*0.5
      Wc(i,j,k) = ( U3(i,j,k) + U3(i,j,k+1) )*0.5
      ENDDO
      ENDDO
      ENDDO
!
      DO j = 1, M2
      DO i = 1, M1
      Wa(i,j,0) = Wa(i,j,M3)   ! set the guard cell centre velocity
      Wa(i,j,N3) = Wa(i,j,1)   ! with periodic B.C. for 2 or 6 point 
      Wb(i,j,0) = Wb(i,j,M3)   ! spatial averaging
      Wb(i,j,N3) = Wb(i,j,1)
      Wc(i,j,0) = Wc(i,j,M3)
      Wc(i,j,N3) = Wc(i,j,1)
      ENDDO
      ENDDO
!
      IF( Model_Dim .EQ. 3 )THEN 
!
      DO k = 1, M3
      DO j = 1, M2
      Wa(0,j,k) = Wa(1,j,k)
      Wa(N1,j,k) = Wa(M1,j,k)
      Wb(0,j,k) = Wb(1,j,k)
      Wb(N1,j,k) = Wb(M1,j,k)
      Wc(0,j,k) = Wc(1,j,k)
      Wc(N1,j,k) = Wc(M1,j,k)
      ENDDO
      ENDDO
!
      DO k = 1, M3
      DO i = 1, M1
      Wa(i,0,k) = Wa(i,1,k)
      Wa(i,N2,k) = Wa(i,M2,k)
      Wb(i,0,k) = Wb(i,1,k)
      Wb(i,N2,k) = Wb(i,M2,k)
      Wc(i,0,k) = Wc(i,1,k)
      Wc(i,N2,k) = Wc(i,M2,k)
      ENDDO
      ENDDO

!
!     Compute the fluctuating velocity from the accumulated mean.
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      Wa(i,j,k) = Wa(i,j,k) - u1_acc(i,j)
      Wb(i,j,k) = Wb(i,j,k) - u2_acc(i,j)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

!   
!     Rdy_T, Rdy_B : Top, Bottom
!     Rdx_L, Rdx_R : Left, Right

!
      IF( Model_Dim .EQ. 3 )THEN
!
      DO j = 1, M2
      Rdy_T(j) = 2.0 / ( dy(j) + dy(j+1) )
      Rdy_B(j) = 2.0 / ( dy(j-1) + dy(j) )
      ENDDO
!
      DO i = 1, M1
      Rdx_L(i) = 2.0 / ( dx(i-1) + dx(i) )
      Rdx_R(i) = 2.0 / ( dx(i) + dx(i+1) )
      ENDDO
!
!     Filtered 2nd-order velocity function contributed by resolved velocity
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1     
      Delta = ( dx(i)*dy(j)*Delta_z )**(OneThird) 
!     
      Wd(i,j,k) = OneSixth * (                                        &
       ( (Delta*Rdx_R(i))**2) * ( (Wa(i,j,k) - Wa(i+1,j,k))**2         &
      + (Wb(i,j,k) - Wb(i+1,j,k))**2 + (Wc(i,j,k) - Wc(i+1,j,k))**2 ) &                             
      + ( (Delta*Rdx_L(i))**2) * ( (Wa(i,j,k) - Wa(i-1,j,k))**2        &
      + (Wb(i,j,k) - Wb(i-1,j,k))**2 + (Wc(i,j,k) - Wc(i-1,j,k))**2 ) &
      + ( (Delta*Rdy_T(j))**2) * ( (Wa(i,j,k) - Wa(i,j+1,k))**2        &
      + (Wb(i,j,k) - Wb(i,j+1,k))**2 + (Wc(i,j,k) - Wc(i,j+1,k))**2 ) &
      + ( (Delta*Rdy_B(j))**2) * ( (Wa(i,j,k) - Wa(i,j-1,k))**2        &
      + (Wb(i,j,k) - Wb(i,j-1,k))**2 + (Wc(i,j,k) - Wc(i,j-1,k))**2 ) &
      + ( (Delta*Rdz)**2) * ( (Wa(i,j,k) - Wa(i,j,k+1))**2             &
      + (Wb(i,j,k) - Wb(i,j,k+1))**2 + (Wc(i,j,k) - Wc(i,j,k+1))**2 ) &
      + ( (Delta*Rdz)**2) * ( (Wa(i,j,k) - Wa(i,j,k-1))**2             &
      + (Wb(i,j,k) - Wb(i,j,k-1))**2 + (Wc(i,j,k) - Wc(i,j,k-1))**2 )  )
      ENDDO
      ENDDO
      ENDDO
!
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
!     Total structure function; correction factor is taken into account for
!     unresolved velocity.
!     Multiplying factor is 0.04044 without correction.
!
      Delta = ( dx(i)*dy(j)*Delta_z )**(OneThird)
      E(i,j,k) = Nu + 0.06338* Delta* SQRT( Wd(i,j,k) )    
      ENDDO
      ENDDO
      ENDDO
!
      ENDIF  ! endif of Model_Dim=3

      IF( Model_Dim .EQ. 1 )THEN
!
!
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      Wd(i,j,k) = 0.5 * (                                      &
      (Wa(i,j,k)-Wa(i,j,k-1))**2 + (Wb(i,j,k)-Wb(i,j,k-1))**2  &               
      + (Wc(i,j,k)-Wc(i,j,k-1))**2  +                          &
      (Wa(i,j,k)-Wa(i,j,k+1))**2 + (Wb(i,j,k)-Wb(i,j,k+1))**2  &
      + (Wc(i,j,k)-Wc(i,j,k+1))**2  )          
      ENDDO
      ENDDO
      ENDDO
!
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      E(i,j,k) = Nu + 0.06338*Delta_z* SQRT( Wd(i,j,k) )
      ENDDO
      ENDDO
      ENDDO

!
      ENDIF  ! endif of Model_Dim=1
!     

      ENDIF  ! endif of structure function model
!-------------------------------------------------------------
!
!     Set eddy viscosity at the corner guard cells.
      DO k = 1, M3
      E(0,0,k)  = E(1,1,k)    ! used in the momentum calculations
      E(N1,0,k) = E(M1,1,k)
      E(0,N2,k) = E(1,M2,k)
      E(N1,N2,k) = E(M1,M2,k)
      ENDDO

!
!     Set eddy viscosity at the guard cells (zero normal gradient B.C.)

      DO k = 1, M3
      DO i = 1, M1
      E(i,0,k)  = E(i,1,k)
      E(i,N2,k) = E(i,M2,k)
      ENDDO
      ENDDO
!

      DO k = 1, M3
      DO j = 1, M2
      E(0,j,k) = E(1,j,k)
      E(N1,j,k) = E(M1,j,k)
      ENDDO
      ENDDO

!
!     Set periodic B.C. for the eddy viscosity.

      DO j = 0, N2
      DO i = 0, N1
      E(i,j,0)  = E(i,j,M3)
      E(i,j,N3) = E(i,j,1)
      ENDDO
      ENDDO


!     
!
      IF( Model .EQ. None )THEN
      DO k = 0, N3
      DO j = 0, N2
      DO i = 0, N1
      E(i,j,k) = Nu
      ENDDO
      ENDDO
      ENDDO
      ENDIF
!
!     Maximum viscous stability number (< 0.3 )
!     Vmax = dt*Max{ (E()*(1/dx**2+1/dy**2+1/dz**2) }
!     

      DO k = 0, N3
      DO j = 0, N2
      DO i = 0, N1
      Wa(i,j,k) = 0.0  ! Otherwise MAXVAL will mistakenly pick up the  
      ENDDO            ! guard cell data stored in Wa()
      ENDDO
      ENDDO

!
max=0.0
      DO k = 1, M3
      DO j = 1, M2
      DO i = 1, M1
      Wa(i,j,k) = Delta_t*( E(i,j,k)*(dx0(i)*dx0(i)+dy0(j)*dy0(j)+  &
                                 Rdz*Rdz)  )
      IF (Wa(i,j,k) .GT. max) THEN
         max = Wa(i,j,k)
          imax=i
          jmax=j
          kmax=k
          ENDIF
      ENDDO
      ENDDO
      ENDDO
      Viscous_no = MAXVAL( Wa )
      Emax = MAXVAL( E )


!      WRITE(*,*) 'Max Viscous Stability No.= ', Viscous_no
!      WRITE(*,*) 'at i=',imax, 'j=',jmax,'k=', kmax
!
      IF( Model .NE. None )THEN
      Ratio = (Emax-Nu)*100.0/Nu
      WRITE(*,*) 'Ratio of Max. Eddy and Molecular Visc. (%) = ', Ratio
      ENDIF
!
      RETURN
      END
 
!-----------------------------------------------------------------------
      SUBROUTINE Temperature
!--------------------------------
!Advance the scalar field
        INCLUDE 'jet3d.inc'
        INTEGER :: i,j,k, idum, jdum, kdum
        REAL ::  sum
        REAL :: TMAX, TMIN
        REAL :: TIME1, TIME2
        
IF(Step_No .EQ. 1) THEN
DO I=0, M1
DO J=0, N2
DO K=0, N3
T1(i,j,k) = T1(1,j,k)
ENDDO
ENDDO
ENDDO
TMAX = 1.0
TMIN = 0.0
ENDIF

!make this less messy for run restarts etc.
TMAX = 1.0
TMIN = 0.0

!b.c's

    DO I=1,M1
    DO K=0,N3
      T1(i,N2,k) = T1(i,M2,k)
      T1(i,0,k) = T1(i,1,k)
    ENDDO
    ENDDO

    DO I=1, M1
    DO J=0, N2
      T1(i,j,N3) = T1(i,j,M3)
      T1(i,j,0) = T1(i,j,1)
    ENDDO
    ENDDO

!Work out Alpha for field

DO K=0,N3
DO J=0, N2
DO I=0, N1
 Wh(i,j,k) = 2.2E-5 + ( E(i,j,k) / 0.3 ) 
ENDDO
ENDDO
ENDDO




!U1 velocity component



DO J=0,N2
DO K=0,N3
Wa(1,j,k) = 0.5* (T1(1,j,k) + T1(0,j,k))
IF (Wa(1,j,k) .GT. TMAX) THEN
Wa(1,j,k) = TMAX
ENDIF
IF (Wa(1,j,k) .LT. TMIN) THEN
Wa(1,j,k) = TMIN
ENDIF
ENDDO
ENDDO

 


DO K=0, N3
DO J=0, N2
DO I=2, N1
!Wa(i,j,k) = 0.5* (T1(i,j,k) + T1(i-1,j,k))
Wa(i,j,k) = 1.5*(T1(i-1,j,k)) - 0.5*(T1(i-2,j,k))
IF (Wa(i,j,k) .GT. TMAX) THEN
Wa(i,j,k) = TMAX
ENDIF
IF (Wa(i,j,k) .LT. TMIN) THEN
Wa(i,j,k) = TMIN
ENDIF
ENDDO
ENDDO
ENDDO



DO K=0, N3
DO J=0, N2
DO I=1, M1
!Wa(i,j,k) = 0.5* (T1(i,j,k) + T1(i-1,j,k))
Wd(i,j,k) = 1.5*(T1(i,j,k)) - 0.5*(T1(i+1,j,k))
IF (Wd(i,j,k) .GT. TMAX) THEN
Wd(i,j,k) = TMAX
ENDIF
IF (Wd(i,j,k) .LT. TMIN) THEN
Wd(i,j,k) = TMIN
ENDIF
ENDDO
ENDDO
ENDDO

 



    DO K=0, N3
    DO J=0, N2
    DO I=0, N1
    IF (U1(i,j,k) .GT. 0.0) THEN
       Wb(i,j,k) = (-U1(i,j,k) * Wa(i,j,k)) + Wh(i,j,k) * dx0(i) * &
             (T1(i,j,k) - T1(i-1,j,k))
    ELSE
       Wb(i,j,k) = (-U1(i,j,k) * Wd(i,j,k)) + Wh(i,j,k) * dx0(i) * &
             (T1(i,j,k) - T1(i-1,j,k))
    ENDIF
    ENDDO
    ENDDO
    ENDDO 





        DO K=0, N3
        DO J=0,N2
        DO I=1, M1
            Wc(i,j,k) = dx0(i) * (Wb(i+1,j,k) - Wb(i,j,k))
        ENDDO
        ENDDO
        ENDDO

     
     

!U2 component

        DO I=1, M1
        DO K=0, N3
          Wa(i,1,k) = 0.5 * (T1(i,1,k) + T1(i,0,k))
            IF (Wa(i,j,k) .GT. TMAX) THEN
            Wa(i,j,k) = TMAX
            ENDIF
            IF (Wa(i,j,k) .LT. TMIN) THEN
            Wa(i,j,k) = TMIN
            ENDIF
        ENDDO
        ENDDO




        DO K=0, N3
        DO J=2, N2
        DO I=1, M1
            !Wa(i,j,k) = 0.5* (T1(i,j,k) + T1(i,j-1,k))
            Wa(i,j,k) = 1.5 * T1(i,j-1,k) - 0.5 *(T1(i,j-2,k))
            IF (Wa(i,j,k) .GT. TMAX) THEN
            Wa(i,j,k) = TMAX
            ENDIF
            IF (Wa(i,j,k) .LT. TMIN) THEN
            Wa(i,j,k) = TMIN
            ENDIF
        ENDDO
        ENDDO
        ENDDO




        DO K=0, N3
        DO J=1, M2
        DO I=1, M1
            Wd(i,j,k) = 1.5*T1(i,j,k) - 0.5*(T1(i,j+1,k))
            IF (Wd(i,j,k) .GT. TMAX) THEN
            Wd(i,j,k) = TMAX
            ENDIF
            IF (Wd(i,j,k) .LT. TMIN) THEN
            Wd(i,j,k) = TMIN
            ENDIF
        ENDDO
        ENDDO
        ENDDO

 



    DO K=0, N3
    DO J=1,N2
    DO I=1, M1
        IF(U2(i,j,k) .GT. 0.0000000) THEN
        Wb(i,j,k) = (-U2(i,j,k) * Wa(i,j,k)) + Wh(i,j,k) * dy0(j) * (T1(i,j,k) - T1(i,j-1,k))
        ELSE
        Wb(i,j,k) = (-U2(i,j,k) * Wd(i,j,k)) + Wh(i,j,k) * dy0(j) * (T1(i,j,k) - T1(i,j-1,k))
        ENDIF
    ENDDO
    ENDDO
    ENDDO

 



    DO K=0, N3
    DO J=0,M2
    DO I=1, M1
        Wc(i,j,k) = Wc(i,j,k) + dy0(j) * (Wb(i,j+1,k) - Wb(i,j,k))
    ENDDO
    ENDDO
    ENDDO 



!U3 component


    DO I=1, M1
    DO J=0, N2
        Wa(i,j,1) = 0.5 * (T1(i,j,1) + T1(i,j,0))
        IF (Wa(i,j,1) .GT. TMAX) THEN
        Wa(i,j,1) = TMAX
        ENDIF
        IF (Wa(i,j,1) .LT. TMIN) THEN
        Wa(i,j,1) = TMIN
    ENDIF
    ENDDO
    ENDDO



    DO K=2, N3
    DO J=0, N2
    DO I=1, M1
        Wa(i,j,k) = 1.5 * T1(i,j,k-1) - 0.5*T1(i,j,k-2) 
        IF (Wa(i,j,k) .GT. TMAX) THEN
        Wa(i,j,k) = TMAX
        ENDIF
        IF (Wa(i,j,k) .LT. TMIN) THEN
        Wa(i,j,k) = TMIN
        ENDIF
    ENDDO
    ENDDO
    ENDDO

    

DO K=1, M3
DO J=0,N2
DO I=1, M1
Wd(i,j,k) = 1.5*T1(i,j,k) - 0.5*T1(i,j,k+1)
IF (Wd(i,j,k) .GT. TMAX) THEN
Wd(i,j,k) = TMAX
ENDIF
IF (Wd(i,j,k) .LT. TMIN) THEN
Wd(i,j,k) = TMIN
ENDIF
ENDDO
ENDDO
ENDDO
    

DO K=1, N3
DO J=0, N2
DO I=1, M1
IF (U3(i,j,k) .GT. 0.000000) THEN
Wb(i,j,k) = (-U3(i,j,k) * Wa(i,j,k)) +  Wh(i,j,k) * Rdz * (T1(i,j,k) - T1(i,j,k-1))
ELSE
Wb(i,j,k) = (-U3(i,j,k) * Wd(i,j,k)) +  Wh(i,j,k) * Rdz * (T1(i,j,k) - T1(i,j,k-1))
ENDIF
ENDDO
ENDDO
ENDDO



DO K=1,M3
DO J=0,N2
DO I=1, M1
Wc(i,j,k) = Wc(i,j,k) + Rdz * (Wb(i,j,k+1) - Wb(i,j,k))
ENDDO
ENDDO
ENDDO


!T update - Adams Bashfoth method

      DO K=1, M3
      DO j = 1, M2 
      DO I=1, M1
        T1(i,j,k) = T1(i,j,k) + Alpha*Wc(i,j,k) + Beta*R1(i,j,k)
        R1(i,j,k) = Wc(i,j,k)
      ENDDO 
      ENDDO 
      ENDDO


IF (Step_No .EQ. 20000) THEN
IF(MOD(Step_No, 100) .EQ. 0) THEN
DO I=1,N1
DO J=0,N2
DO K=0, N3
IF (T1(i,j,k) .GT. TMAX) THEN
T1(i,j,k) = TMAX
ENDIF
IF (T1(i,j,k) .LT. TMIN) THEN
T1(i,j,k) = TMIN
ENDIF
ENDDO
ENDDO
ENDDO
ENDIF
ENDIF

      RETURN

      END SUBROUTINE Temperature
!------------------------------------------------------
      SUBROUTINE Euler_PDF
!--------------------------------
!Calculate scalar fields using Eulerian PDF method
        INCLUDE 'jet3d.inc'
		INCLUDE 'mpif.h'
        INTEGER :: i,j,k, idum, jdum, kdum, root,yplate
        REAL ::  sum, maxval

		REAL :: CO1(0:(N1+1)*(N2+1)*(N3+1)),CO2(0:(N1+1)*(N2+1)*(N3+1)),CO3(0:(N1+1)*(N2+1)*(N3+1))
		REAL :: CO4(0:(N1+1)*(N2+1)*(N3+1)),CO5(0:(N1+1)*(N2+1)*(N3+1))
		REAL :: CA1(0:(N1+1)*(N2+1)*(N3+1)),CA2(0:(N1+1)*(N2+1)*(N3+1)),CA3(0:(N1+1)*(N2+1)*(N3+1))
		REAL :: CA4(0:(N1+1)*(N2+1)*(N3+1)),CA5(0:(N1+1)*(N2+1)*(N3+1))
        REAL :: ftmp(0:N1,0:N2,0:N3),randno(1:9),rand


IF(nprocs .NE. nfields) THEN
write(*,*) 'Number of cpus used does not match number of stochastic fields.'
Abort = .TRUE.
ENDIF


CALL MPI_BCAST(Step_No, 1, MPI_INTEGER,0, MPI_COMM_WORLD,ierr)

!Initialisation - needs to go somewhere else too WAM 27/01/2013
IF(Step_No .EQ. 1) THEN


   do k=0, n3
    do j=0, n2
     if(passstream .eq. 1) then
     if(j .ge. jsplit) then
   pass1(1,j,k) = 1.0
   pasf1(1,j,k) = pass1(1,j,k)
   pasf1(0,j,k) =pasf1(1,j,k)
     else
   pass1(1,j,k) = 0.0
   pasf1(1,j,k) = pass1(1,j,k)
   pasf1(0,j,k) =pasf1(1,j,k)
    endif !j
    else if(passstream .eq. 2) then
   if(j .ge. jsplit) then
   pass1(1,j,k) = 0.0
   pasf1(1,j,k) = pass1(1,j,k)
   pasf1(0,j,k) =pasf1(1,j,k)
   else
   pass1(1,j,k) = 1.0
   pasf1(1,j,k) = pass1(1,j,k)
   pasf1(0,j,k) =pasf1(1,j,k)
   endif !j
   endif ! passstream
enddo
   enddo

	DO I = 2, N1
	DO J=1, M2
	DO K=1, M3
	pass1(i,j,k) = pass1(1,j,k)
	ENDDO
	ENDDO
	ENDDO

	DO I=1, N1
	DO J=1, M2
	DO K=1, M3
	pasf1(i,j,k) = pasf1(1,j,k)
	ENDDO
	ENDDO
	ENDDO
	
	DO I=0, N1
	DO J=0, N2
	pasf1(i,j,0) = pasf1(i,j,M3)
	pasf1(i,j,N3) = pasf1(i,j,1)
	ENDDO
	ENDDO

	DO I=0,N1
	DO K=0,N3
	pasf1(i,0,k) = pasf1(i,1,k)
	pasf1(i,N2,k) = pasf1(i,M2,k)
	ENDDO
	ENDDO

   
    do j = 0, N2
    do k = 0, N3
     if(NOstream .eq. 1) then
      if( j .ge. jsplit) then
      rNO(1,j,k) = NOin
      rNO(0,j,k) = NOin
      rO3(1,j,k) = 0.0
      rO3(0,j,k) = 0.0
      NOf(1,j,k) = NOin
      NOf(0,j,k) = NOin
      O3f(1,j,k) = 0.0
      O3f(0,j,k) = 0.0
     else
      rNO(1,j,k) = 0.0
      rNO(0,j,k) = 0.0
      rO3(1,j,k) = O3in
      rO3(0,j,k) = O3in
      NOf(1,j,k) = 0.0
      NOf(0,j,k) = 0.0
      O3f(1,j,k) = O3in
      O3f(0,j,k) = O3in
     endif
     else ! nostream
     if( j .ge. jsplit) then
     rNO(1,j,k) = 0.0
     rNO(0,j,k) = 0.0
     rO3(1,j,k) = O3in
     rO3(0,j,k) = O3in
     NOf(1,j,k) = 0.0
     NOf(0,j,k) = 0.0
     O3f(1,j,k) = O3in
     O3f(0,j,k) = O3in
     else
     rNO(1,j,k) = NOin
     rNO(0,j,k) = NOin
     rO3(1,j,k) = 0.0
     rO3(0,j,k) = 0.0
     NOf(1,j,k) = NOin
     NOf(0,j,k) = NOin
     O3f(1,j,k) = 0.0
     O3f(0,j,k) = 0.0
     endif ! j
     endif !nostream
    enddo
    enddo

    do i = 2, n1
    do J=0, n2
    do k=0, n3
     rNO(i,j,k) = rNO(1,j,k)
     rO3(i,j,k) = rO3(1,j,k)
     NOf(i,j,k) = NOf(1,j,k)
     O3f(i,j,k) = O3f(1,j,k)
    enddo
    enddo
    enddo

ENDIF


IF(myrank .EQ. 0) THEN
!pack up required variables into single dimension arrays
    DO I=0, N1
    DO K=0, N3
    DO J=0, N2
        CO1( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = U1(i,j,k)
        CO2( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = U2(i,j,k)
        CO3( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = U3(i,j,k)
        CO4( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = E(i,j,k)
    ENDDO
    ENDDO
    ENDDO
ENDIF

!Broadcast them
CALL MPI_BCAST(CO1, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(CO2, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(CO3, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(CO4, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)

!and then unpack
IF(myrank .GT. 0) THEN
    DO I=0, N1
    DO K=0, N3
    DO J=0, N2
        U1(i,j,k) = CO1( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
        U2(i,j,k) = CO2( i*((N2+1)*(N3+1)) + k*(N2+1) +j) 
        U3(i,j,k) = CO3( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
        E(i,j,k) =  CO4( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
    ENDDO
    ENDDO
    ENDDO
ENDIF

DO idum = 0, (N1+1)*(N2+2)*(N3+1)
    CO1(idum) = 0.0
    CO2(idum) = 0.0
    CO3(idum) = 0.0
    CO4(idum) = 0.0
ENDDO

!IF(myrank .EQ. 0) THEN
!WRITE(*,*) 'Values broadcast'
!ENDIF


!Call random numbers for the Weiner process here WAM 22/04/13
do i=1,9
call random_number(rand)
rand = 2.0*( rand - 0.5 )
if(rand .gt. 0.0) randno(i) = 1.0
if(rand .le. 0.0) randno(i) = -1.0
enddo

!do i=1, 9
!write(910+myrank,*) step_no, randno(i)
!enddo


!Possible to split micromixing term out of calc_scal
! to a new subroutine at the end of the process?



! This set of calls is for the passive scalar
!First stage
CALL calc_scal(pasf1,pass1,Wc,0,0.0, Wi,Wi,0)
call update_scal(Wa,pasf1,Wc,1,Wc)
CALL setbc_scal(Wa,0)
!Second stage
CALL calc_scal(Wa,pass1,Wd,0,0.0, Wi,Wi,0)
call update_scal(Wb,pasf1,Wd,2, Wa)
call setbc_scal(Wb,0)
!Third stage
CALL CALC_SCAL(Wb,pass1,We,0,0.0,Wi,Wi,0)
call update_scal(Wg,pasf1,We,3, Wb)
CALL setbc_SCAL(Wg,0)
!call weiner process here
call weiner_proc(ftmp, Wg,randno(1),randno(2),randno(3))

DO I=0, N1
DO K=0, N3
DO J=0, N2
CO1( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = ftmp(i,j,k)
pasf1(i,j,k) = ftmp(i,j,k)
ENDDO
ENDDO
ENDDO
call setbc_scal(pasf1,0)

! This is the set of calls for the reacting species
!  NO + O3 -> NO2 + O2
!
!Stage 1 - NO
call calc_scal(NOf, rNO,We,1,rrate1,O3f,Wi,0)
call update_scal(Wa,NOf,We,1,Wc)
call setbc_scal(Wa,1)

!Stage 1 - O3
call calc_scal(O3f, rO3,We,1,rrate2,NOf,Wi,0)
call update_scal(Wb,O3f,We,1,Wc)
call setbc_scal(Wb,2)

!Stage 1 - NO2
call calc_scal(NO2f, rNO2,We,1,rrate3,NOf, O3f,1)
call update_scal(Wc,NO2f,We,1,Wc)
call setbc_scal(Wc,3)

!Stage 1 completed

!Stage 2 - NO
call calc_scal(Wa,rNO,We,1,rrate1,Wb,Wi,0)
call update_scal(Wa,NOf,We,2,Wa)
call setbc_scal(Wa,1)

!Stage 2 - O3
call calc_scal(Wb,rO3,We,1,rrate2,Wa,Wi,0)
call update_scal(Wb,O3f,We,2,Wb)
call setbc_scal(Wb,2)

!Stage 2 - NO2
call calc_scal(Wc,rNO2,We,1,rrate3,Wa,Wb,1)
call update_scal(Wc,NO2f,We,2,Wc)
call setbc_scal(Wc,3)

!Stage 2 completed

!Stage 3 - NO
call calc_scal(Wa,rNO,We,1,rrate1,Wb,Wi,0)
call update_scal(Wa,NOf,We,3,Wa)
call setbc_scal(Wa,1)
call weiner_proc(ftmp,Wa,randno(1),randno(2),randno(3))

DO I=0, N1
DO K=0, N3
DO J=0, N2
CO2( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = ftmp(i,j,k)
NOf(i,j,k) = ftmp(i,j,k)
ENDDO
ENDDO
ENDDO
call setbc_scal(NOf,1)


!Stage 3 - O3
call calc_scal(Wb,rO3,We,1,rrate2,Wa,Wi,0)
call update_scal(Wb,O3f,We,3,Wb)
call setbc_scal(Wb,2)
call weiner_proc(ftmp,Wb,randno(1),randno(2),randno(3))

DO I=0, N1
DO K=0, N3
DO J=0, N2
CO3( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = ftmp(i,j,k)
O3f(i,j,k) = ftmp(i,j,k)
ENDDO
ENDDO
ENDDO
call setbc_scal(O3f,2)

!Stage 3 - NO2
call calc_scal(Wc,rNO2,We,1,rrate3,Wa, Wb,1)
call update_scal(Wc,NO2f,We,3,Wc)
call setbc_scal(Wc,3)
call weiner_proc(ftmp,Wc,randno(1),randno(2),randno(3))

DO I=0, N1
DO K=0, N3
DO J=0, N2
CO4( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = ftmp(i,j,k)
NO2f(i,j,k) = ftmp(i,j,k)
ENDDO
ENDDO
ENDDO
call setbc_scal(NO2f,3)


if(mod(step_no, 1000) .eq. 0) then
WRITE(400+myrank,*) 'VARIABLES = "x" "Y" "pass1" "NO" "O3" "NO2"'
WRITE(400+myrank,*) 'ZONE I=', M2,', J=', M1
DO I=1, M1
DO J=1,M2
WRITE(400+myrank,*) x(i), y(j), pasf1(i,j,1),NOf(i,j,1),O3f(i,j,1),NO2f(i,j,1)
ENDDO
ENDDO
endif
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!Change this code to reduce onto the master, and then broadcast to other processes
!this should save communication time
CALL MPI_REDUCE(CO1, CA1, (N1+1)*(N2+1)*(N3+1), MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(CO2, CA2, (N1+1)*(N2+1)*(N3+1),MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(CO3, CA3, (N1+1)*(N2+1)*(N3+1),MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(CO4, CA4, (N1+1)*(N2+1)*(N3+1),MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(CO5, CA5, (N1+1)*(N2+1)*(N3+1),MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)


    if(myrank .eq. 0) then
    DO I=0, N1
    DO K=0, N3
    DO J=0, N2
CA1( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = CA1( i*((N2+1)*(N3+1)) + k*(N2+1) +j)/real(nprocs)
CA2( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = CA2( i*((N2+1)*(N3+1)) + k*(N2+1) +j)/real(nprocs)
CA3( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = CA3( i*((N2+1)*(N3+1)) + k*(N2+1) +j)/real(nprocs)
CA4( i*((N2+1)*(N3+1)) + k*(N2+1) +j) = CA4( i*((N2+1)*(N3+1)) + k*(N2+1) +j)/real(nprocs)
        pass1(i,j,k) = CA1( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
        rNO(i,j,k) = CA2( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
        rO3(i,j,k) = CA3( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
        rNO2(i,j,k) = CA4( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
    ENDDO
    ENDDO
    ENDDO
    endif

CALL MPI_BCAST(CA1, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(CA2, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(CA3, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(CA4, (N1+1)*(N2+1)*(N3+1), MPI_REAL,0, MPI_COMM_WORLD,ierr)

!Broadcast scalar fields to other processes
IF(myrank .GT. 0) THEN
DO I=0, N1
DO K=0, N3
DO J=0, N2
pass1(i,j,k) = CA1( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
rNO(i,j,k) = CA2( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
rO3(i,j,k) = CA3( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
rNO2(i,j,k) =  CA4( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
rO2(i,j,k) = CA5( i*((N2+1)*(N3+1)) + k*(N2+1) +j)
ENDDO
ENDDO
ENDDO
ENDIF

DO idum = 0, (N1+1)*(N2+2)*(N3+1)
CO1(idum) = 0.0
CO2(idum) = 0.0
CO3(idum) = 0.0
CO4(idum) = 0.0
CO5(idum) = 0.0
CA1(idum) = 0.0
CA2(idum) = 0.0
CA3(idum) = 0.0
CA4(idum) = 0.0
CA5(idum) = 0.0
ENDDO


!set conditions
DO K=0, N3
DO J=0, N2
DO I=0, N1
  if (pass1(i,j,k) .LT. 0.0) pass1(i,j,k) = 0.0
  if (pass1(i,j,k) .GT. 1.0) pass1(i,j,k) = 1.0
if (rNO(i,j,k) .LT. 0.0) rNO(i,j,k) = 0.0
if (rNO(i,j,k) .GT. NOin) rNO(i,j,k) = NOin
if (rO3(i,j,k) .LT. 0.0) rO3(i,j,k) = 0.0
if (rO3(i,j,k) .GT. O3in) rO3(i,j,k) = O3in
if (rNO2(i,j,k) .LT. 0.0) rNO2(i,j,k) = 0.0
if (rO2(i,j,k) .LT. 0.0) rO2(i,j,k) = 0.0
ENDDO
ENDDO
ENDDO


CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
      RETURN

      END SUBROUTINE Euler_PDF

!------------------------------------------------
       subroutine update_scal(output,xi,operat,stage,vars)
!This subroutine updates the scalar in steps
!stage corresponds to the stage of the runge-kutta scheme
    include 'jet3d.inc'

integer :: i,j,k, stage
real :: operat(0:n1,0:n2,0:n3), xi(0:n1,0:n2,0:n3), output(0:n1,0:n2,0:n3)
real :: vars(0:n1,0:n2,0:n3)
real :: oneoverthree
real :: twooverthree


oneoverthree = 1.0/3.0
twooverthree = 2.0/3.0

do i=1,m1
do j=1,m2
do k=1,m3
! put wiener process in here
if(stage .eq. 1) then
output(i,j,k) = xi(i,j,k) + Delta_t*operat(i,j,k) 
else if(stage .eq. 2) then
output(i,j,k) = 0.75*xi(i,j,k) + 0.25*vars(i,j,k) + 0.25*Delta_t*operat(i,j,k)
else if(stage .eq. 3) then
output(i,j,k) =  oneoverthree*xi(i,j,k) + twooverthree*vars(i,j,k) + twooverthree*Delta_t*operat(i,j,k)
endif
enddo
enddo
enddo

return

end subroutine
!---------------------------------------------------------------
       subroutine setbc_SCAL(work,var)
! sets bounday conditions for the stochastic fields
! var indicates which field is being used
! 0 = passive scalar, 1 = NO, 2 = O3, 3 = NO2, 4 = O2
include 'jet3d.inc'

integer :: i,j,k,var
real :: work(0:n1,0:n2,0:n3)

DO I=0, N1
DO J=0, N2
work(i,j,0) = work(i,j,M3)
work(i,j,N3) = work(i,j,1)
ENDDO
ENDDO

DO I=0,N1
DO K=0,N3
work(i,0,k) = work(i,1,k)
work(i,N2,k) = work(i,M2,k)
ENDDO
ENDDO

do j=0,n2
do k=0, n3
!outflow bc
work(n1,j,k) = work(m1,j,k)
!inflow bcs for scalar fields
! 0 is passive scalar
if(var .eq. 0) then
if(passstream .eq. 1) then
if (j .lt. jsplit) then
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
else
!work(1,j,k) = 1.0
work(0,j,k) = 1.0
endif
else if(passstream .eq. 2) then
if (j .lt. jsplit) then
!work(1,j,k) = 1.0
work(0,j,k) = 1.0
else
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
endif
endif ! passstream
endif !var

! 1 is NO
if(var .eq. 1) then
if(NOStream .eq. 1) then
if(j .lt. jsplit) then
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
else
!work(1,j,k) = NOin
work(0,j,k) = NOin
endif !j
else !nostream
if(j .lt. jsplit) then
!work(1,j,k) = NOin
work(0,j,k) = NOin
else
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
endif !j
endif !Nostream
endif !var

!O3
if(var .eq. 2) then
if(NOStream .eq. 1) then
if(j .lt. jsplit) then
!work(1,j,k) = O3in
work(0,j,k) = O3in
else
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
endif !j
else !no stream
if(j .lt. jsplit) then
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
else
!work(1,j,k) = O3in
work(0,j,k) = O3in
endif !j
endif !Nostream
endif !var

!products
if(var .eq. 3 .or. var .eq. 4) then
!work(1,j,k) = 0.0
work(0,j,k) = 0.0
endif

enddo
enddo

do i=0,n1
do j=0,n2
do k=0,n3
if(work(i,j,k) .lt. 0.0) work(i,j,k) = 0.0
if(var .eq. 1) then
if(work(i,j,k) .gt. NOin) work(i,j,k) = NOin
endif
if(var .eq. 2) then
if(work(i,j,k) .gt. O3in) work(i,j,k) = O3in
endif
enddo
enddo
enddo


return

end subroutine
!-----------------------------------------------------------
      subroutine weiner_proc(in, xi,ran1,ran2,ran3)

include 'jet3d.inc'

integer :: i,j,k
real :: in(0:n1,0:n2,0:n3), xi(0:n1,0:n2,0:n3)
real :: rand, fef,fwf,ran1,ran2,ran3

! apply Weiner process to computation of new stochastic field
! need a dichotomic vector, and each random number has to be differet
! on each cpu

!weiner term is (2*Gamma')**0.5 * (d xi)/(dx) dW
! Gamma' = mu/sigma + mu_sgs/sigma_sgs
! sigma = schmidt number, sigma_sgs = sgs schmidt number, 0.3 - 0.7

! weiner process is different for each field and each spatial co-ordinate
! but independent of spatial location x

! i.e. random number for each field and each term in (d xi / dx), applied 
! over all of the domain

call random_number(rand)
rand = 2.0*( rand - 0.5 )
if(rand .gt. 0.0) rand = 1.0
if(rand .le. 0.0) rand = -1.0

!x-direction
do i = 1, m1
do j = 1, m2
do k = 1, m3
  fef = 0.5*(xi(i,j,k) + xi(i+1,j,k))
  fwf = 0.5*(xi(i,j,k) + xi(i-1,j,k))
in(i,j,k) = xi(i,j,k) +dx0(i)*(fef-fwf)*ran1*sqrt(Delta_t) * &
            sqrt(2*( (Nu/ sclam) + (E(i,j,k) - Nu) / scturb))
enddo
enddo
enddo

call random_number(rand)
rand = 2.0*( rand - 0.5 )
if(rand .gt. 0.0) rand = 1.0
if(rand .le. 0.0) rand = -1.0

!y-direction
do i = 1, m1
do j = 1, m2
do k = 1, m3
fef = 0.5*(xi(i,j,k) + xi(i,j+1,k))
fwf = 0.5*(xi(i,j,k) + xi(i,j-1,k))
in(i,j,k) = in(i,j,k) + dy0(j)*(fef-fwf)*ran2*sqrt(Delta_t) * &
            sqrt(2*( (Nu/ sclam) + (E(i,j,k) - Nu) / scturb))
enddo
enddo
enddo

call random_number(rand)
rand = 2.0*( rand - 0.5 )
if(rand .gt. 0.0) rand = 1.0
if(rand .le. 0.0) rand = -1.0

!z-direction
do i = 1, m1
do j = 1, m2
do k = 1, m3
fef = 0.5*(xi(i,j,k) + xi(i,j,k+1))
fwf = 0.5*(xi(i,j,k) + xi(i,j,k-1))
in(i,j,k) = in(i,j,k) + Rdz*(fef-fwf)*ran3*sqrt(Delta_t)  * &
            sqrt(2*( (Nu/ sclam) + (E(i,j,k) - Nu) / scturb))
enddo
enddo
enddo


return

end subroutine
!--------------------------------------------------------------
SUBROUTINE CALC_SCAL(scal,ph,Waa,react,rate,reactant1,reactant2,product)

        INCLUDE 'jet3d.inc'
		INCLUDE 'mpif.h'
        REAL :: small_num,rate
		REAL :: scal(0:N1,0:N2,0:N3) 
        REAL :: Work, ph(0:N1,0:N2,0:N3), Waa(0:n1,0:n2,0:n3)
        real :: reactant1(0:N1,0:N2,0:N3), reactant2(0:N1,0:N2,0:N3)
		INTEGER :: i,j,k, react, product
		REAL :: dfdxeef, dfdxef, &
		         dfdxwf, dfdxwwf, &
				 dfdyssf, dfdysf, &
				 dfdynf, dfdynnf, &
				 dfdzllf, dfdzlf, &
				 dfdzrf, dfdzrrf
		 REAL :: ratio
		 REAL :: lime, limw,      &
		         lims, limn,      &
				 liml, limr
		REAL :: fef, fwf,         &
		        fnf, fsf,         &
				flf, frf
		REAL:: gwf,gef,           &
		       gsf,gnf,           &
			   glf,grf
		REAL :: galwf,galef,     &
		       galsf, galnf,      &
			   gallf,galrf
		REAL :: twf, tef,         &
			   tsf, tnf,          &
			   tlf, trf
		REAL :: scalfwf, scalfef, &
		      scalfsf, scalfnf,   &
			  scalflf,scalfrf

small_num = 1.0E-30

DO K=1, M3
DO J=1, M2
DO I=1, M1
!ef,wf,nf,sf,lf,rf
   dfdxef = ( scal(i+1,j,k) - scal(i,j,k) ) / dxv(i)
   dfdxwf = ( scal(i,j,k) - scal(i-1,j,k) ) / dxv(i-1)
   dfdynf = ( scal(i,j+1,k) - scal(i,j,k) ) / dyu(j)
   dfdysf = ( scal(i,j,k) - scal(i,j-1,k) ) / dyu(j-1)
   dfdzrf = ( scal(i,j,k+1) - scal(i,j,k) ) / ( z(k+1) - z(k) )
   dfdzlf = ( scal(i,j,k) - scal(i,j,k-1) ) / ( z(k) - z(k-1) )
!eef,wwf,nnf,ssf,llf,rrf
   dfdxeef = ( scal(i+2,j,k) - scal(i+1,j,k) ) / dxv(i+1)
   dfdxwwf = ( scal(i-1,j,k) - scal(i-2,j,k) ) / dxv(i-2)
   dfdynnf = ( scal(i,j+2,k) - scal(i,j+1,k) ) / dyu(j+1)
   dfdyssf = ( scal(i,j-1,k) - scal(i,j-2,k) ) / dyu(j-2)
   dfdzrrf = ( scal(i,j,k+2) - scal(i,j,k+1) ) / ( z(k+2) - z(k+1) )
   dfdzllf = ( scal(i,j,k-1) - scal(i,j,k-2) ) / ( z(k-1) - z(k-2) )
IF( I .EQ. 1) THEN
   dfdxwwf = dfdxwf
ENDIF
IF( I .EQ. M1) THEN
   dfdxeef = dfdxef
ENDIF
IF( J .EQ. 1) THEN
   dfdyssf = dfdysf
ENDIF
IF( J .EQ. M2) THEN
   dfdynnf = dfdynf
ENDIF
IF( K .EQ. 1) THEN
   dfdzllf = dfdzlf
ENDIF
IF( K .EQ. M3) THEN
   dfdzrrf = dfdzrf
ENDIF

!gradient ratio calculator

			IF(U1(i,j,k) .GE. 0.0) THEN
				ratio = (dfdxwwf + small_num) / (dfdxwf + small_num)
			ELSE
				ratio = (dfdxef + small_num) / (dfdxwf + small_num)
			ENDIF

                        limw = amax1(0.0,amin1(2*ratio,1.0))

			IF(U1(i+1,j,k) .GT. 0.0) THEN
			   ratio = (dfdxwf + small_num) / (dfdxef + small_num)
			ELSE
			   ratio = (dfdxeef + small_num) / (dfdxef + small_num)
			ENDIF

			lime = amax1(0.0,amin1(2*ratio,1.0))

			IF(U2(i,j,k) .GE. 0.0) THEN
			   ratio = (dfdyssf + small_num) / (dfdysf + small_num)
			ELSE
			   ratio = (dfdynf + small_num) / (dfdysf + small_num)
			ENDIF

                        lims = amax1(0.0,amin1(2*ratio,1.0))
			
			IF(U2(i,j+1,k) .GT. 0.0) THEN
			   ratio = (dfdysf + small_num) / (dfdynf + small_num)
			ELSE
			   ratio = (dfdynnf + small_num) / (dfdynf + small_num)
			ENDIF

                        limn = amax1(0.0,amin1(2*ratio,1.0))

			IF(U3(i,j,k) .GE. 0.0) THEN
			   ratio = (dfdzllf + small_num) / (dfdzlf + small_num)
			ELSE
			   ratio = (dfdzrf + small_num) / (dfdzlf + small_num)
			ENDIF
           
                          liml = amax1(0.0,amin1(2*ratio,1.0))

			IF(U3(i,j,k+1) .GT. 0.0) THEN
			   ratio = (dfdzlf + small_num) / (dfdzrf + small_num)
			ELSE
			   ratio = (dfdzrrf + small_num) / (dfdzrf + small_num)
			ENDIF
		       
                        limr = amax1(0.0,amin1(2*ratio,1.0))


!TVD Flux limiter:  Branley and Jones


		IF(U1(i,j,k) .GT. 0.0) THEN
		   fwf = scal(i-1,j,k) + 0.5*limw*(scal(i,j,k) - scal(i-1,j,k))
		ELSE IF(U1(i,j,k) .LT. 0.0) THEN
		   fwf = scal(i,j,k) + 0.5*limw*(scal(i-1,j,k) - scal(i,j,k))
		ELSE
		   fwf = 0.5*(scal(i-1,j,k) + scal(i,j,k))
		ENDIF

		IF(U1(i+1,j,k) .GT. 0.0) THEN
		   fef = scal(i,j,k) + 0.5*lime*(scal(i+1,j,k) - scal(i,j,k))
		ELSE IF(U1(i+1,j,k) .LT. 0.0) THEN
		   fef = scal(i+1,j,k) + 0.5*lime*(scal(i,j,k) - scal(i+1,j,k))
		ELSE
		   fef = 0.5*(scal(i,j,k) + scal(i+1,j,k))
		ENDIF

		IF(U2(i,j,k) .GT. 0.0) THEN
		   fsf = scal(i,j-1,k) + 0.5*lims*(scal(i,j,k) - scal(i,j-1,k))
		ELSE IF(U2(i,j,k) .LT. 0.0) THEN
		   fsf = scal(i,j,k) + 0.5*lims*(scal(i,j-1,k) - scal(i,j,k))
		ELSE
		   fsf = 0.5*(scal(i,j,k) + scal(i,j-1,k))
		ENDIF

		IF(U2(i,j+1,k) .GT. 0.0) THEN
		  fnf = scal(i,j,k) + 0.5*limn*(scal(i,j+1,k) - scal(i,j,k))
		ELSE IF(U2(i,j+1,k) .LT. 0.0) THEN
		  fnf = scal(i,j+1,k) + 0.5*limn*(scal(i,j,k) - scal(i,j+1,k))
		ELSE
		  fnf = 0.5*(scal(i,j,k) + scal(i,j+1,k))
		ENDIF

		IF(U3(i,j,k) .GT. 0.0) THEN
		   flf = scal(i,j,k-1) + 0.5*liml*(scal(i,j,k) - scal(i,j,k-1))
		ELSE IF(U3(i,j,k) .LT. 0.0) THEN
		   flf = scal(i,j,k) + 0.5*liml*(scal(i,j,k-1) - scal(i,j,k))
		ELSE
		   flf = 0.5*(scal(i,j,k) + scal(i,j,k-1))
		ENDIF

		IF(U3(i,j,k+1) .GT. 0.0) THEN
		   frf = scal(i,j,k) + 0.5*limr*(scal(i,j,k+1) - scal(i,j,k))
		ELSE IF(U3(i,j,k+1) .LT. 0.0) THEN
		   frf = scal(i,j,k+1) + 0.5*limr*(scal(i,j,k) - scal(i,j,k+1))
		ELSE
		   frf = 0.5*(scal(i,j,k) + scal(i,j,k+1))
		ENDIF
		
		IF(I .EQ. 1) THEN
		fwf = 0.5*(scal(i,j,k) + scal(i-1,j,k))
		fef = 0.5*(scal(i,j,k) + scal(i+1,j,k))
		ENDIF
		
		IF(I .EQ. M1) THEN
		fwf = 0.5*(scal(i,j,k) + scal(i-1,j,k))
		fef = 0.5*(scal(i,j,k) + scal(i+1,j,k))
		ENDIF


!all fluxes computed - get convective fluxes in cells  = ui*psi

  gwf = U1(i,j,k)*fwf*areawf(i,j,k)
  gef = U1(i+1,j,k)*fef*areaef(i,j,k)
  gsf = U2(i,j,k)*fsf*areasf(i,j,k)
  gnf = U2(i,j+1,k)*fnf*areanf(i,j,k)
  glf = U3(i,j,k)*flf*arealf(i,j,k)
  grf = U3(i,j,k+1)*frf*arealf(i,j,k)

 Work = -(gef - gwf) - (gnf - gsf) - (grf - glf)

!diffusion fluxes terms


  galwf = (0.5 * (  ((E(i,j,k)-Nu)/Scturb) + ((E(i-1,j,k)-Nu)/Scturb)) + Nu/Sclam)
  galef = (0.5 * (  ((E(i+1,j,k)-Nu)/Scturb) + ((E(i,j,k)-Nu)/Scturb)) + Nu/Sclam)
  galsf = (0.5 * (  ((E(i,j,k)-Nu)/Scturb) + ((E(i,j-1,k)-Nu)/Scturb)) + Nu/Sclam)
  galnf = (0.5 * (  ((E(i,j+1,k)-Nu)/Scturb) + ((E(i,j,k)-Nu)/Scturb)) + Nu/Sclam)
  gallf = (0.5 * (  ((E(i,j,k)-Nu)/Scturb) + ((E(i,j,k-1)-Nu)/Scturb)) + Nu/Sclam)
  galrf = (0.5 * (  ((E(i,j,k+1)-Nu)/Scturb) + ((E(i,j,k)-Nu)/Scturb)) + Nu/Sclam)


  twf = galwf * areawf(i,j,k) * dfdxwf
  tef = galef * areaef(i,j,k) * dfdxef
  tsf = galsf * areasf(i,j,k) * dfdysf
  tnf = galnf * areanf(i,j,k) * dfdynf
  tlf = gallf * arealf(i,j,k) * dfdzlf
  trf = galrf * arearf(i,j,k) * dfdzrf

  Work = Work + (tef - twf) + (tnf - tsf) + (trf - tlf)

  Work = Work / (dx(i)*dy(j)*dz(k))

! now need to add in mixing and reaction terms
! move this into a separate fractional step?

    galwf = (E(i,j,k)-Nu)/Scturb + (Nu/Sclam)
    Waa(i,j,k) = Work - (0.5 * (galwf / ( (dx(i)*dy(j)*dz(k))**(2.0/3.0))) * (scal(i,j,k) - ph(i,j,k)))

 ! reaction terms 
    if(react .eq. 1 .and. product .eq. 1) then
    Waa(i,j,k) = Waa(i,j,k) + rate*reactant1(i,j,k)*reactant2(i,j,k)
    else
    Waa(i,j,k) = Waa(i,j,k) + rate*scal(i,j,k)*reactant1(i,j,k)
    endif

ENDDO
ENDDO
ENDDO


RETURN

END SUBROUTINE CALC_SCAL
!--------------------------------end of subroutine--------------------------------------

       SUBROUTINE Momentum
!      -------------------
!      Advance the momentum vector field, provisionally
       INCLUDE 'jet3d.inc'
 
       INTEGER i, j, k
 
!       IF LES then Nu -> E(i,j,k) & accumulate s^2 in Wa(i,j,k)
 
!      U1 equation


       DO i = 0, M1 
       DO j = 0, M2 
       DO k = 0, N3
!        Strain & Stress11 = Wc
         Wc(i,j,k) = (U1(i+1,j,k)-U1(i,j,k))/dx(i)
         Wc(i,j,k) = 2.0 * E(i,j,k) * Wc(i,j,k)                     &
                   - 0.25*(U1(i+1,j,k)+U1(i,j,k))**2
       ENDDO 
       ENDDO
       ENDDO



       DO i = 1, N1 
       DO j = 0, N2 
       DO k = 0, N3
!        Wa = part sum of Force 1
         Wa(i,j,k) = dx1(i) * ( Wc(i,j,k) - Wc(i-1,j,k) )    
       ENDDO 
       ENDDO
       ENDDO



       DO i = 1, N1 
       DO j = 1, N2 
       DO k = 0, N3
!       Strain & Stress12: Wb ~= s21
        Wb(i,j,k) = (U2(i,j,k) - U2(i-1,j,k))*dx1(i)                &
                  + (U1(i,j,k) - U1(i,j-1,k)) * dy1(j)
        Wd(i,j,k) = E(i,j,k) * Wb(i,j,k) - 0.25 *                   &
        (U2(i,j,k) + U2(i-1,j,k)) * (U1(i,j,k) + U1(i,j-1,k))
      ENDDO 
      ENDDO
      ENDDO



      DO i = 0, N1
      DO j = 0, M2 
      DO k = 0, N3
        Wa(i,j,k) = Wa(i,j,k) + dy0(j) * ( Wd(i,j+1,k)-Wd(i,j,k) )  
      ENDDO 
      ENDDO 
      ENDDO



      DO i = 1, N1 
      DO j = 0, N2 
      DO k = 1, N3
!       St = Strain & Stress 13
        Wd(i,j,k) = Rdz*(U1(i,j,k) - U1(i,j,k-1))                    &
                  + dx1(i)*(U3(i,j,k) - U3(i-1,j,k))
        Wd(i,j,k) = E(i,j,k) * Wd(i,j,k)                             &
                  - 0.25 * (U1(i,j,k) + U1(i,j,k-1))                 &
                         * (U3(i,j,k) + U3(i-1,j,k))
      ENDDO 
      ENDDO 
      ENDDO



      DO i = 0, N1 
      DO j = 0, N2 
      DO k = 0, M3
         Wa(i,j,k) = Wa(i,j,k) + Rdz * ( Wd(i,j,k+1) - Wd(i,j,k))
      ENDDO 
      ENDDO
       ENDDO



      DO i = 0, M1 
      DO j = 0, M2 
      DO k = 0, N3
!       Wd = Strain and Stress 22
        Wd(i,j,k) = dy0(j) * (U2(i,j+1,k) - U2(i,j,k))
        Wd(i,j,k) = 2.0 * E(i,j,k) * Wd(i,j,k)                      &
                    - 0.25 * (U2(i,j+1,k) + U2(i,j,k))**2
      ENDDO 
      ENDDO 
      ENDDO


      DO i = 1, N1 
      DO j = 0, N2 
      DO k = 0, N3
        Wa(i,j,k) = Wa(i,j,k) + Force_x
      ENDDO 
      ENDDO 
      ENDDO

!     U2 equation ***************
 

      DO i = 1, N1
      DO j = 1, N2 
      DO k = 0, N3
!       Stress21 using s12 saved in Wb
        Wb(i,j,k) = E(i,j,k) * Wb(i,j,k)                               &
                  - 0.25 * ( U2(i,j,k) + U2(i-1,j,k) )                 &
                  * ( U1(i,j,k) + U1(i,j-1,k) )
      ENDDO 
      ENDDO 
      ENDDO



      DO i = 0, M1 
      DO j = 1, N2 
      DO k = 0, N3
        Wb(i,j,k) = dy1(j)*(Wd(i,j,k) -Wd(i,j-1,k))                    &
                  + (Wb(i+1,j,k)-Wb(i,j,k) ) / dx(i) 
      ENDDO 
      ENDDO
      ENDDO
!     Note that use of Wb here is non-recursive as forward-looking.


      DO i = 0, N1 
      DO j = 1, N2 
      DO k = 1, N3
!       Strain and Stress23 = Wc; strain saved in St
        Wd(i,j,k) = Rdz*(U2(i,j,k) - U2(i,j,k-1)) + dy1(j)*            &
                    (U3(i,j,k) - U3(i,j-1,k))
        Wc(i,j,k) = E(i,j,k) * Wd(i,j,k) - 0.25 *                      &
        (U2(i,j,k)+U2(i,j,k-1) ) * (U3(i,j,k)+U3(i,j-1,k))
      ENDDO 
      ENDDO
      ENDDO



     
      DO i = 0, N1 
      DO j = 0, N2 
      DO k = 0, M3
        Wb(i,j,k) = Wb(i,j,k) + Rdz*(Wc(i,j,k+1)-Wc(i,j,k)) + Force_y
      ENDDO 
      ENDDO
      ENDDO


!     U3 equation ***************
 

      DO i = 1, N1 
      DO j = 0, N2 
      DO k = 1, N3
!       Strain & Stress31
        Wc(i,j,k) = Rdz*(U1(i,j,k) - U1(i,j,k-1)) +                    &
                  dx1(i)*(U3(i,j,k) - U3(i-1,j,k))
        Wc(i,j,k) = E(i,j,k) * Wc(i,j,k)                               &
                  - 0.25 * (U1(i,j,k) + U1(i,j,k-1))                   &
                         * (U3(i,j,k) + U3(i-1,j,k))                     
      ENDDO 
      ENDDO 
      ENDDO


      DO i = 0, M1
      DO j = 0, N2 
      DO k = 0, N3
        Wc(i,j,k) = (Wc(i+1,j,k) - Wc(i,j,k))/dx(i)
!       not recursive as forward looking
      ENDDO 
      ENDDO
      ENDDO


      DO i = 0, N1
      DO j = 1, N2 
      DO k = 1, N3
!       Stress32
        Wd(i,j,k) = E(i,j,k) * Wd(i,j,k)                               &
                  - 0.25 * (U2(i,j,k)+U2(i,j,k-1))                     &
                  * (U3(i,j,k) + U3(i,j-1,k))
      ENDDO
      ENDDO 
      ENDDO



      DO i = 0, N1 
      DO j = 0, M2 
      DO k = 0, N3
        Wc(i,j,k) = Wc(i,j,k) + dy0(j)*(Wd(i,j+1,k) - Wd(i,j,k))
      ENDDO 
      ENDDO
      ENDDO



      DO i = 0, N1 
      DO j = 0, N2 
      DO k = 0, M3
!       Strain & Stress33
        Wd(i,j,k) = 2.0 * E(i,j,k) * Rdz * (U3(i,j,k+1) - U3(i,j,k))   &
                  - 0.25 * (U3(i,j,k+1) + U3(i,j,k))**2
      ENDDO 
      ENDDO 
      ENDDO



      DO i = 0, N1
      DO j = 0, N2 
      DO k = 1, N3
        Wc(i,j,k) = Wc(i,j,k) + Rdz*(Wd(i,j,k) - Wd(i,j,k-1))          &
                  + Force_z
      ENDDO 
      ENDDO 
      ENDDO

 
!     U1 update

      DO i = 1, M1
      DO j = 1, M2 
      DO k = 1, M3
        U1(i,j,k) = U1(i,j,k) + Alpha*Wa(i,j,k) + Beta*H1(i,j,k)
        H1(i,j,k) = Wa(i,j,k)
      ENDDO 
      ENDDO 
      ENDDO


!     U2 update

      DO i = 1, M1
      DO j = 1, M2 
      DO k = 1, M3
        U2(i,j,k) = U2(i,j,k) + Alpha*Wb(i,j,k) + Beta*H2(i,j,k)
        H2(i,j,k) = Wb(i,j,k)
      ENDDO 
      ENDDO
      ENDDO


!     U3 update

      DO i = 1, M1
      DO j = 1, M2 
      DO k = 1, M3
        U3(i,j,k) = U3(i,j,k) + Alpha*Wc(i,j,k) + Beta*H3(i,j,k)
        H3(i,j,k) = Wc(i,j,k)
      ENDDO 
      ENDDO
      ENDDO

 
      RETURN
      END

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE Pressure
!     -------------------
!     Fourier/MGADZSLR hybrid Poisson solver.
      INCLUDE 'jet3d.inc'
 
      INTEGER i, j, k
      REAL sum,P_red


!     DIVERGENCE = source term
      sum = 0.0
      DO i = 1, M1 
      DO j = 1, M2 
      DO k = 1, M3
       Wa(i,j,k) = twoverdt * ( (U1(i+1,j,k) - U1(i,j,k)) * dx0(i)     &
                 + dy0(j) * (U2(i,j+1,k) - U2(i,j,k))                  &
                 + Rdz * (U3(i,j,k+1) - U3(i,j,k)) )
       sum = sum + Wa(i,j,k)
      ENDDO 
      ENDDO 
      ENDDO

!     Ui are the fluxes, not velocities if not Cartesian
 
      CALL LUD3(Nor1,Sou1,Est1,Wst1,A1,Wb,Wc,NxA,NyA,1)
!     do every time to set Wb, Wc work arrays = B1, C1 coefficients

      CALL fft_Forward(Wa)
      CALL fft_Forward(P1)

      IF (cells .EQ. tall) THEN
        DO k = 1, M3
          CALL NS_Zebra(k)
        ENDDO
      ENDIF
      IF (cells .EQ. wide) THEN
        DO k = 1, M3
          CALL EW_Zebra(k)
        ENDDO
      ENDIF
      IF (cells .EQ. mixt) THEN
        DO k = 1, M3
          CALL AD_Zebra(k)
        ENDDO
      ENDIF
 
      CALL fft_Inverse(P1) 

!     Setting Pressure Boundary Condition, Generell
      CALL Pressure_bc

!     Correct the Velocity Field
      DO i = 1, M1 
      DO j = 1, M2 
      DO k = 1, M3
        U1(i,j,k) = U1(i,j,k) - dtover2 * dx1(i)                       &
                  * (P1(i,j,k) - P1(i-1,j,k))
      ENDDO 
      ENDDO 
      ENDDO

      DO i = 1, M1 
      DO j = 1, M2 
      DO k = 1, M3
        U2(i,j,k) = U2(i,j,k) - dtover2 * dy1(j)                       &
                  * (P1(i,j,k) - P1(i,j-1,k))
      ENDDO 
      ENDDO 
      ENDDO

      DO i = 1, M1 
      DO j = 1, M2 
      DO k = 1, M3
        U3(i,j,k) = U3(i,j,k)                                          &
                  - dtover2 * Rdz * (P1(i,j,k) - P1(i,j,k-1))
      ENDDO 
      ENDDO 
      ENDDO
 
!     Reduce the Pressure Field by P1(1,1,1)
      P_red=P1(1,1,1)
      DO i=0, N1
      DO j=0, N2
      DO k=0, N3
         P1(i,j,k)=P1(i,j,k)-P_red     !#changed
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE NS_Zebra(K)
!     ----------------------
!     2D MultiGrid ZSLR. NS solutions, EW relaxation.
      INCLUDE 'jet3d.inc'
      INTEGER  K,Iter1,Iter2,Iter3,Iter4,Imax1,Imax2,Imax3,Imax4
      PARAMETER (Imax1=25,Imax2=200,Imax3=150,Imax4=25)

      max_residual = 1.0
      threshold = 1.0
      Iter1=0
      Iter2=0
 
      CALL Solve_NS(P1,Wa,Nor1,Sou1,Est1,Wst1,Wb,G1,NxA,NyA,K)
      IF (max_residual .LT. epsilon) GOTO 50

 10   CONTINUE
         Iter1 = Iter1 + 1
         CALL Solve_NS(P1,Wa,Nor1,Sou1,Est1,Wst1,Wb,G1,NxA,NyA,K)
      IF (Iter1 .LE. Imax1 .AND. max_residual .GE. epsilon             &
          .AND. ksqu(O3+K) .GE. 0.16) GOTO 10
      IF (max_residual .LT. epsilon) GOTO 50

!     Use multigrid if not converged:
 20   CONTINUE

         Iter2 = Iter2 + 1
         threshold = max_residual/10.0
         CALL Inject(P1,Wa,P2,R2,Nor1,Sou1,Est1,Wst1,A1,NxB,NyB,K)
         CALL Solve_NS(P2,R2,Nor2,Sou2,Est2,Wst2,B2,G2,NxB,NyB,K)
         CALL Inject(P2,R2,P3,R3,Nor2,Sou2,Est2,Wst2,A2,NxC,NyC,K)

         Iter3 = 0
 30      CONTINUE
            Iter3 = Iter3 + 1
            CALL Solve_NS(P3,R3,Nor3,Sou3,Est3,Wst3,B3,G3,NxC,NyC,K)
         IF (Iter3 .LE. Imax3 .AND. max_residual .GE. threshold        &
             .AND. max_residual .GE. roundoff) GOTO 30

         CALL Interp(P2,P3,NxC,NyC,K)
         CALL Solve_NS(P2,R2,Nor2,Sou2,Est2,Wst2,B2,G2,NxB,NyB,K)
         CALL Solve_NS(P2,R2,Nor2,Sou2,Est2,Wst2,B2,G2,NxB,NyB,K)
         CALL Interp(P1,P2,NxB,NyB,K)
         CALL Solve_NS(P1,Wa,Nor1,Sou1,Est1,Wst1,Wb,G1,NxA,NyA,K)
         IF (max_residual .LT. epsilon) GOTO 50

         Iter4 = 0
 40      CONTINUE
            Iter4 = Iter4 + 1
            CALL Solve_NS(P1,Wa,Nor1,Sou1,Est1,Wst1,Wb,G1,NxA,NyA,K)
         IF (Iter4 .LE. Imax4 .AND. max_residual .LE. 1.1*epsilon      &
             .AND. max_residual .GE. epsilon) GOTO 40

      IF ( max_residual .GE. epsilon .AND. Iter2 .LE. Imax2 ) GOTO 20
 
 50   CONTINUE

      IF (Iter2 .GT. Imax2 .AND. Debug .GT. 0) THEN
            WRITE (*,*)'Iter1-4=',Iter1,Iter2,Iter3,Iter4    
            WRITE (*,*)'Max. Res.=',max_residual 
            WRITE (*,*)'->WARNING: non convergent in NS_Zebra!'
      ENDIF

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE EW_Zebra(K)
!     ----------------------
!     2D MultiGrid ZSLR. EW solutions, NS relaxation.
      INCLUDE 'jet3d.inc'
      INTEGER K,Iter1,Iter2,Iter3,Iter4,Imax1,Imax2,Imax3,Imax4
      PARAMETER (Imax1=25,Imax2=200,Imax3=150,Imax4=25) 

      max_residual = 1.0
      threshold = 1.0
      Iter1=0
      Iter2=0

      CALL Solve_EW(P1,Wa,Nor1,Sou1,Est1,Wst1,Wc,G1,NxA,NyA,K)
      IF (max_residual .LT. epsilon) GOTO 50

 10   CONTINUE
         Iter1 = Iter1 + 1
         CALL Solve_EW(P1,Wa,Nor1,Sou1,Est1,Wst1,Wc,G1,NxA,NyA,K)
      IF (Iter1.LE.Imax1.AND.max_residual.GE.epsilon.AND.ksqu(O3+K)    &
          .GE.0.16) GOTO 10
      IF (max_residual .LT. epsilon) GOTO 50

!     Use multigrid if not converged: 
 20   CONTINUE

         Iter2 = Iter2 + 1
         threshold = max_residual/10.0
         CALL Inject(P1,Wa,P2,R2,Nor1,Sou1,Est1,Wst1,A1,NxB,NyB,K)
         CALL Solve_EW(P2,R2,Nor2,Sou2,Est2,Wst2,C2,G2,NxB,NyB,K)
         CALL Inject(P2,R2,P3,R3,Nor2,Sou2,Est2,Wst2,A2,NxC,NyC,K)

         Iter3=0
 30      CONTINUE
            Iter3 = Iter3 + 1
            CALL Solve_EW(P3,R3,Nor3,Sou3,Est3,Wst3,C3,G3,NxC,NyC,K)
         IF (max_residual.GE.threshold.AND.max_residual.GE.roundoff    &
             .AND.Iter3.LE.Imax3) GOTO 30

         CALL Interp(P2,P3,NxC,NyC,K)
         CALL Solve_EW(P2,R2,Nor2,Sou2,Est2,Wst2,C2,G2,NxB,NyB,K)
         CALL Solve_EW(P2,R2,Nor2,Sou2,Est2,Wst2,C2,G2,NxB,NyB,K)
         CALL Interp(P1,P2,NxB,NyB,K)
         CALL Solve_EW(P1,Wa,Nor1,Sou1,Est1,Wst1,Wc,G1,NxA,NyA,K)
         IF (max_residual .LT. epsilon) GOTO 50

         Iter4=0
 40      CONTINUE
            Iter4 = Iter4 + 1
            CALL Solve_EW(P1,Wa,Nor1,Sou1,Est1,Wst1,Wc,G1,NxA,NyA,K)
         IF (max_residual.LE.1.1*epsilon.AND.max_residual.GE.epsilon   &
             .AND.Iter4.LE.Imax4) GOTO 40

      IF (max_residual .GE. epsilon .AND. Iter2 .LE. Imax2) GOTO 20

 50   CONTINUE

      IF (Iter2 .GT. Imax2 .AND. Debug .GT. 0) THEN
            WRITE (*,*)'Iter1-4=',Iter1,Iter2,Iter3,Iter4    
            WRITE (*,*)'Max. Res.=',max_residual 
            WRITE (*,*)'->WARNING: non convergent in EW_Zebra!'
      ENDIF

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE AD_Zebra(K)
!     ----------------------
!     2D MultiGrid ZSLR. Alternating direction.
      INCLUDE 'jet3d.inc'
      INTEGER K,Iter1,Iter2,Iter3,Iter4,Imax1,Imax2,Imax3,Imax4  
      PARAMETER (Imax1=25,Imax2=200,Imax3=150,Imax4=25)

      max_residual = 1.0
      threshold = 1.0
      Iter1=0
      Iter2=0

 10   CONTINUE
         Iter1 = Iter1 + 1
         CALL Solve_NS(P1,Wa,Nor1,Sou1,Est1,Wst1,Wb,G1,NxA,NyA,K)
         CALL Solve_EW(P1,Wa,Nor1,Sou1,Est1,Wst1,Wc,G1,NxA,NyA,K)
      IF (Iter1 .LE. Imax1 .AND. max_residual .GE. epsilon             &
          .AND. ksqu(O3+K) .GE. 0.16 ) GOTO 10
      IF (max_residual .LT. epsilon) GOTO 50

!     Use multigrid if not converged: 
 20   CONTINUE

         Iter2 = Iter2 + 1
         threshold = max_residual/10.0
         CALL Inject(P1,Wa,P2,R2,Nor1,Sou1,Est1,Wst1,A1,NxB,NyB,K)
         CALL Solve_NS(P2,R2,Nor2,Sou2,Est2,Wst2,B2,G2,NxB,NyB,K)
         CALL Solve_EW(P2,R2,Nor2,Sou2,Est2,Wst2,C2,G2,NxB,NyB,K)
         CALL Inject(P2,R2,P3,R3,Nor2,Sou2,Est2,Wst2,A2,NxC,NyC,K)

         Iter3 = 0
 30      CONTINUE
            Iter3 = Iter3 + 1
            CALL Solve_NS(P3,R3,Nor3,Sou3,Est3,Wst3,B3,G3,NxC,NyC,K)
            CALL Solve_EW(P3,R3,Nor3,Sou3,Est3,Wst3,C3,G3,NxC,NyC,K)
         IF (Iter3 .LE. Imax3 .AND. max_residual.GE.threshold          &
             .AND. max_residual .GE. roundoff) GOTO 30

         CALL Interp(P2,P3,NxC,NyC,K)
         CALL Solve_NS(P2,R2,Nor2,Sou2,Est2,Wst2,B2,G2,NxB,NyB,K)
         CALL Solve_EW(P2,R2,Nor2,Sou2,Est2,Wst2,C2,G2,NxB,NyB,K)
         CALL Interp(P1,P2,NxB,NyB,K)

         Iter4 = 0
 40      CONTINUE
            Iter4 = Iter4 + 1
            CALL Solve_NS(P1,Wa,Nor1,Sou1,Est1,Wst1,Wb,G1,NxA,NyA,K)
            CALL Solve_EW(P1,Wa,Nor1,Sou1,Est1,Wst1,Wc,G1,NxA,NyA,K)
         IF (Iter4 .LE. Imax4 .AND. max_residual .LE. 1.1*epsilon      &
             .AND. max_residual .GE. epsilon) GOTO 40

      IF (max_residual .GE. epsilon .AND. Iter2 .LE. Imax2) GOTO 20

 50   CONTINUE

      IF (Iter2 .GT. Imax2 .AND. Debug .GT. 0) THEN
            WRITE (*,*)'Iter1-4  =',Iter1,Iter2,Iter3,Iter4    
            WRITE (*,*)'Max. Res.=',max_residual
            WRITE (*,*)'->WARNING: non convergent in AD_Zebra!'
      ENDIF

      RETURN
      END 

!----------------------------------------------------------------------
      SUBROUTINE Solve_NS(P,R,Nor,Sou,Est,Wst,B,G,Nx,Ny,K)
!     ----------------------------------------------------
!     Forward and backward substitution in I; zebra striping in J. K fixed.
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'
      INTEGER M, I, J, K, Nx, Ny, Imax, Jmax

      REAL P(0:Nx,0:Ny,0:N3), R(0:Nx,0:Ny,0:N3)
      REAL Nor(0:Nx,0:Ny), Sou(0:Nx,0:Ny)
      REAL Est(0:Nx,0:Ny), Wst(0:Nx,0:Ny)
      REAL B(0:Nx,0:Ny,0:N3), G(0:Nx,0:Ny)
      REAL res

      Imax = Nx-1
      Jmax = Ny-1
      max_residual = 0.0

      DO M = 1, 2
         DO J = M, Jmax, 2
          DO I = 1, Imax
          G(I,J) = R(I,J,K) -P(I,J-1,K)*Est(I,J) -P(I,J+1,K)*Wst(I,J)
          Wd(I,J,K) = P(I,J,K)
          ENDDO
        ENDDO

          DO J = M, Jmax, 2
          DO I = 1, Imax
            G(I,J) = (G(I,J) - G(I-1,J) * Nor(I,J)) * B(I,J,K)
          ENDDO 
          ENDDO

          DO J = M, Jmax, 2
          DO I = Imax, 1,-1
            P(I,J,K) = G(I,J) -Sou(I,J)*B(I,J,K)*P(I+1,J,K)
          ENDDO
        ENDDO !J
      ENDDO !M

      DO J = 1, Jmax
      DO I = 1, Imax
        res = ABS( P(I,J,K) - Wd(I,J,K) )
        IF (res .GT. max_residual) max_residual = res
      ENDDO
      ENDDO

      CALL Press_Planes( P, Nx, Ny, N3, K) !not needed

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE Solve_EW(P,R,Nor,Sou,Est,Wst,C,G,Nx,Ny,K)
!     ----------------------------------------------------
!     Forward and backward substitution in J; zebra striping in I. K fixed.
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'
      INTEGER M, I, J, K, Nx, Ny, Imax, Jmax

      REAL P(0:Nx,0:Ny,0:N3), R(0:Nx,0:Ny,0:N3)
      REAL Nor(0:Nx,0:Ny), Sou(0:Nx,0:Ny)
      REAL Est(0:Nx,0:Ny), Wst(0:Nx,0:Ny)
      REAL C(0:Nx,0:Ny,0:N3), G(0:Nx,0:Ny)
      REAL temp, res

      Imax = Nx-1
      Jmax = Ny-1
      max_residual = 0.0
      DO I = 1, Imax
        G(I,0) = 0.0
      ENDDO
  DO M = 1, 2

      DO I = M, Imax, 2
        DO J = 1, Jmax
        G(I,J) = ( R(I,J,K) -P(I-1,J,K)*Nor(I,J) -P(I+1,J,K)*Sou(I,J)  &
                  -G(I,J-1)*Est(I,J) ) * C(I,J,K)
        ENDDO
      ENDDO

      DO I = M, Imax, 2
        DO J = Jmax, 1,-1
         temp = P(I,J,K)
         P(I,J,K) = G(I,J) - Wst(I,J)*C(I,J,K)*P(I,J+1,K)
         res = ABS(P(I,J,K)-temp)
         IF (res .GT. max_residual) max_residual = res
        ENDDO
      ENDDO
     
  ENDDO
      
      CALL Press_Planes( P, Nx, Ny, N3, K)   !not needed

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE Inject(PF,RF,P,R,NF,SF,EF,WF,AF,Nx,Ny,K)
!     ---------------------------------------------------
      INCLUDE 'jet3d.inc'
      INTEGER I, J, K, II, JJ, Nx, Ny, IN, JN
      REAL PF(0:Nx*2-1,0:Ny*2-1,0:N3), RF(0:Nx*2-1,0:Ny*2-1,0:N3) 
      REAL P(0:Nx,0:Ny,0:N3), R(0:Nx,0:Ny,0:N3)
      REAL NF(0:Nx*2-1,0:Ny*2-1),  SF(0:Nx*2-1,0:Ny*2-1) 
      REAL EF(0:Nx*2-1,0:Ny*2-1),  WF(0:Nx*2-1,0:Ny*2-1)
      REAL AF(0:Nx*2-1,0:Ny*2-1)

      IN = Nx-1
      JN = Ny-1
      DO J = 1, JN
      DO I = 1, IN
        P(I,J,K) = 0.0
        R(I,J,K) = 0.0
      ENDDO
      ENDDO


      DO JJ = 1, JN
      DO II = 1, IN
        DO J = JJ*2-1, JJ*2
        DO I = II*2-1, II*2
         R(II,JJ,K) = R(II,JJ,K) + RF(I,J,K)                           &
                    - PF(I,J-1,K)*EF(I,J)-PF(I,J+1,K)*WF(I,J)          &
                    - PF(I+1,J,K)*SF(I,J)-PF(I-1,J,K)*NF(I,J)          &
                    - PF(I,J,K)*(AF(I,J)-ksqu(K))
        ENDDO
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE Interp(Pfine,Pcoarse,Nx,Ny,K)
!     ----------------------------------------
      INCLUDE 'jet3d.inc'
      INTEGER I, J, K, II, JJ, Nx, Ny, IN, JN, Imax, Jmax
      REAL ad, adi, adj, adp, adm
      REAL Pfine(0:Nx*2-1,0:Ny*2-1,0:N3), Pcoarse(0:Nx,0:Ny,0:N3)

      IN = Nx-1
      JN = Ny-1
      DO I= 1, IN
        II=2*I
      DO J= 1, JN
        JJ=2*J
        ad  = Pcoarse(I,J,K)
        adi = (Pcoarse(I+1,J,K)-Pcoarse(I-1,J,K))/4
        adj = (Pcoarse(I,J+1,K)-Pcoarse(I,J-1,K))/4
        adp = ad+adi
        adm = ad-adi
        Pfine(II-1,JJ-1,K) = Pfine(II-1,JJ-1,K)+adm-adj
        Pfine(II-1,JJ,  K) = Pfine(II-1,JJ,  K)+adm+adj
        Pfine(II,  JJ-1,K) = Pfine(II,  JJ-1,K)+adp-adj
        Pfine(II,  JJ,  K) = Pfine(II,  JJ,  K)+adp+adj
      ENDDO
      ENDDO
      Imax = 2*IN
      Jmax = 2*JN

      CALL Press_Planes( Pfine, Nx*2-1, Ny*2-1, N3, K )  !not needed

      RETURN
      END

!----------------------------------------------------------------------
      
      SUBROUTINE MG_table   
!     -------------------
      INCLUDE 'jet3d.inc'
         
         CALL Coeff(Nor1,Sou1,Est1,Wst1,A1,NxA,NyA)
         CALL Coeff(Nor2,Sou2,Est2,Wst2,A2,NxB,NyB)
         CALL Coeff(Nor3,Sou3,Est3,Wst3,A3,NxC,NyC)
!         call for Nor1 etc done in Pressure, every step.         
         CALL LUD3(Nor2,Sou2,Est2,Wst2,A2,B2,C2,NxB,NyB,4)
         CALL LUD3(Nor3,Sou3,Est3,Wst3,A3,B3,C3,NxC,NyC,16)

         RETURN
       END

!----------------------------------------------------------------------
 
      SUBROUTINE Coeff(Nor,Sou,Est,Wst,A,Nx,Ny)
!     -----------------------------------------
!     Mesh coefficients for pressure solver only
      INCLUDE 'jet3d.inc'
      INTEGER Nx, Ny, Imax, Jmax, i, j, idiv, ii, jdiv, jj
      REAL Nor(0:Nx,0:Ny), Sou(0:Nx,0:Ny) 
      REAL Est(0:Nx,0:Ny), Wst(0:Nx,0:Ny)
      REAL A(0:Nx,0:Ny)
      REAL AR
      REAL deltax(0:N1), barx(N1)
      REAL deltay(0:N2), bary(N2)
 
      Imax = Nx -1
      Jmax = Ny -1
      ARmin = 10000000
      ARmax = -ARmin
 
      idiv = M1 / Imax
      DO i = 1, Imax
        deltax(i) = 0.0
        DO ii = 1, idiv
          deltax(i) = deltax(i) + dx((i-1)*idiv+ii)
        ENDDO
      ENDDO
      deltax(0)      = deltax(1)       ! PVM dangerous!
      deltax(Imax+1) = deltax(Imax)
      DO j = 1, Jmax
        jj = j * idiv
        DO i = 1, Imax+1
          barx(i) = 0.5*(deltax(i-1) + deltax(i))
        ENDDO
        DO i = 1, Imax
          Nor(i,j) = idiv**2 / barx(i)   / deltax(i)
          Sou(i,j) = idiv**2 / barx(i+1) / deltax(i)
        ENDDO
        IF(idiv .GT. 1) THEN
        DO i = 1, Imax+1
          barx(i) = 0.5*(deltax(i-1) + deltax(i))
        ENDDO
        DO i = 1, Imax
          Nor(i,j) = idiv**2 / barx(i)   / deltax(i)
          Sou(i,j) = idiv**2 / barx(i+1) / deltax(i)
        ENDDO
        END IF
 
        Nor(1,j) = 0.0
        Sou(Imax,j) = 0.0
      ENDDO
 
      jdiv = M2 / Jmax
      DO j = 1, Jmax
        deltay(j) = 0.0
        DO jj = 1, jdiv
          deltay(j) = deltay(j) + dy((j-1)*jdiv+jj)
        ENDDO
      ENDDO
      deltay(0) = deltay(1)
      deltay(Jmax+1) = deltay(Jmax)
      DO j = 1, Jmax+1
        bary(j) = 0.5 * (deltay(j-1) + deltay(j))
      ENDDO
      DO i = 1, Imax
        ii = i*jdiv
        DO j = 1, Jmax
        Est(i,j) = jdiv**2   /bary(j)   /deltay(j)
        Wst(i,j) = jdiv**2   /bary(j+1) /deltay(j)
        IF(jdiv .GT. 2) THEN
          Est(i,j) = jdiv**2 /bary(j)   /deltay(j)
          Wst(i,j) = jdiv**2 /bary(j+1) /deltay(j)
        END IF
        ENDDO
        Est(i,1) = 0.0
        Wst(i,Jmax) = 0.0
      ENDDO
 
      DO i = 2, Imax
        DO j = 2, Jmax
            AR = Est(i,j)/Nor(i,j)
            IF (AR .LT. ARmin) ARmin=AR
            IF (AR .GT. ARmax) ARmax=AR
        ENDDO
      ENDDO
      DO j = 1, Jmax
        DO i = 1, Imax
          A(i,j) = -Est(i,j)-Wst(i,j)-Nor(i,j)-Sou(i,j)
        ENDDO
      ENDDO
      IF (Imax .LT. M1)  RETURN
      IF (ARmin .LE. 1.0 .AND. ARmax .GT. 1.0) cells = mixt
      IF (ARmax .LE. 1.0 ) cells = tall
      IF (ARmin .GT. 1.0 ) cells = wide

      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE LUD3( Nor, Sou, Est, Wst, A, B, C, Nx, Ny, mult )
!     ------------------------------------------------------------
!     Does LU decomposition of Ny*N3 triadiagonal matrices

      INCLUDE 'jet3d.inc'
      INTEGER I, J, K, Nx, Ny, Imax, Jmax, mult
      REAL Ksq, bb, cc, f, nb, sb, eb, ww, oldsf, oldwf
      REAL Nor(0:Nx,0:Ny), Sou(0:Nx,0:Ny) 
      REAL Est(0:Nx,0:Ny), Wst(0:Nx,0:Ny)
      REAL A(0:Nx,0:Ny), B(0:Nx,0:Ny,0:N3), C(0:Nx,0:Ny,0:N3)
 
      Imax = Nx - 1
      Jmax = Ny - 1
      DO K = 1, M3
!PVM                                Big problems!!!
        Ksq = ksqu(O3+K) * mult
        DO J = 1, Jmax
          B(1,J,K) = 1.0/(A(1,J)-Ksq)
          oldsf = Sou(1,J)*B(1,J,K)                !PVM wrong!
          DO I = 2, Imax
            bb  = 1.0 / (A(I,J)-Ksq)
            nb = Nor(I,J)*bb
            sb = Sou(I,J)*bb
            f  = 1.0 / (1.0 - nb*oldsf)
            B(I,J,K) = bb*f
            oldsf = sb*f
          ENDDO
        ENDDO
 
        DO I = 1, Imax
          C(I,1,K) = 1.0/(A(I,1)-Ksq)
          oldwf = Wst(I,1)*C(I,1,K)      !PVM wrong!
          DO J = 2, Jmax
            cc = 1.0 / (A(I,J)-Ksq)
            eb = Est(I,J) * cc
            ww = Wst(I,J) * cc
            f  = 1.0 / (1.0 - eb*oldwf)
            C(I,J,K) = cc * f
            oldwf = ww * f
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
 
       SUBROUTINE Velocity_bc (set_out)
!      --------------------------------
!      Sets guard planes on velocity
 
         INCLUDE 'jet3d.inc'
         INTEGER bound, oppos
         LOGICAL set_out

!        NOTE: the number of the boundaries with cyclic b.c. must
!        be the largest, & outflow boundary number 2!! 
!PVM     additional routines Heard for PVM interprocessor boundaries
         DO bound = 1, Boundaries
           oppos = bound - (-1)**bound ! goes 2, 1, 4, 3, 6, 5 
           CALL Guard( U1, bound, 1 )
           CALL Heard( U1, oppos, 1 ) 
           CALL Guard( U2, bound, 2 )
           CALL Heard( U2, oppos, 2 )
           CALL Guard( U3, bound, 3 )
           CALL Heard( U3, oppos, 3 )
         END DO

!        Outflow Boundary must be Number 2!
         IF ((Type(2,1) .EQ. Convective) .AND. set_out) CALL Outflow_bc


         RETURN
       END

!----------------------------------------------------------------    

       SUBROUTINE Pressure_bc
!      ----------------------

       INCLUDE 'jet3d.inc'
       INTEGER i,j,k

!      Setting Pressure Boundary Conditions for Each k-Plane
       DO k=1, M3
          CALL Press_Planes (P1,N1,N2,N3,k)
       ENDDO

!      Setting Periodic Pressure BC for the Total Domain in z-Dir.
       DO i=0,N1                              !orig. 1,M1
       DO j=0,N2                              !orig. 1,M2 
            P1(i,j,0)=P1(i,j,M3)
            P1(i,j,N3)=P1(i,j,1)
       ENDDO
       ENDDO

       RETURN
       END

!----------------------------------------------------------------------
 
       SUBROUTINE Guard( A, bo, which )
!      --------------------------------
 
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         INTEGER bo, which, direct
         REAL A(0:N1,0:N2,0:N3)
         INTEGER Jump(3), Hop(3), i, j, k, k1, k2, ir, jr

           k1 = Klo(bo)
           k2 = Khi(bo)

         DO k = 1,3
           Hop(k) = Orient(bo,k) ! step offsets
           IF (Hop(k) .NE. 0) direct = k ! CPVM
           Jump(k) = Period(k) * Hop(k)  ! cyclic vector offsets
         ENDDO
 
!          Given value at boundary
           IF ( Type(bo, which) .EQ. Surface ) THEN
             DO i = Ilo(bo), Ihi(bo)  !orig. Ilo(bo),Ihi(bo)
             DO j = Jlo(bo), Jhi(bo)
             DO k = k1, k2
               A( i-Hop(1), j-Hop(2), k-Hop(3) )                       &
               = bc(bo,which) * 2.0  - A( i, j, k)   
             ENDDO
             ENDDO
             ENDDO
           ENDIF
 
!          ---- inflow boundary (not for normal velocity)
           IF ( Type(bo, which) .EQ. Fixed ) THEN
             DO i = Ilo(bo), Ihi(bo) !orig. Ilo(bo),Ihi(bo)   
             DO j = Jlo(bo), Jhi(bo)
             DO k = k1, k2
               A( i-Hop(1), j-Hop(2), k-Hop(3)) = bc(bo,which)
             ENDDO
             ENDDO
             ENDDO
           ENDIF

!          --- normal velocity = value, on surface
           IF ( Type(bo, which) .EQ. Normal) THEN
             IF ( Orient(bo,which) .EQ. 1 ) THEN ! surface
               Hop(1) = 0
               Hop(2) = 0
               Hop(3) = 0
             ENDIF
             DO i = Ilo(bo), Ihi(bo)  !orig. Ilo(bo),Ihi(bo)
             DO j = Jlo(bo), Jhi(bo)
             DO k = k1, k2
               A( i-Hop(1), j-Hop(2), k-Hop(3) ) = bc(bo,which)
             ENDDO
             ENDDO
             ENDDO
           ENDIF

!          Given gradient*delta-n at boundary
           IF ( Type(bo, which) .EQ. Gradient ) THEN
             DO i = Ilo(bo), Ihi(bo)  !orig. Ilo(bo),Ihi(bo)
             DO j = Jlo(bo), Jhi(bo)
             DO k = k1, k2
               A( i-Hop(1), j-Hop(2), k-Hop(3) )                       &
               = gr(bo,which) + A( i, j, k)   
             ENDDO
             ENDDO
             ENDDO
           ENDIF
 
!          Given sliced inflow
           IF ( Type(bo, which) .EQ. Inflow ) CALL Inflow_bc(bo,which)
 
!          Cyclic b.c. (should not be needed for PVM)
           IF ( Type(bo, which) .EQ. Cyclic ) THEN
!             IF ( bo .EQ. 1 .AND. which .EQ. 1 ) Hop(1) = 0   ! surface
!             IF ( bo .EQ. 6.AND.which .EQ. 3 ) Jump(3) = Jump(3) + 1   !surface
             DO i = Ilo(bo), Ihi(bo)          !orig. Ilo(bo),Ihi(bo)   
             DO j = Jlo(bo), Jhi(bo)          !orig. Jlo(bo),Jhi(bo) 
             DO k = k1, k2
               A( i-Hop(1), j-Hop(2), k-Hop(3) )                       &
               =  A( i+Jump(1), j+Jump(2), k+Jump(3))   
             ENDDO
             ENDDO
             ENDDO
           ENDIF
          
!          Interprocessor talk b.c : sending data
           IF ( Type(bo, which) .EQ. PVM_talk ) THEN
!            Send A( i, j, k ) to next proc
              PRINT *, 'PVM called for ', bo, which
             CALL pvmfinitsend( PVMDEFAULT, info )
             msgtype = which*10 + direct*Hop(direct) ! unique code

           IF ( direct .EQ. 1) THEN
             j  = Jlo(bo)
             i  = Ilo(bo)
             jr = Jhi(bo) - Jlo(bo) + 1
             DO k = k1, k2
               CALL pvmfpack( REAL8, A(i,j,k), jr, N2+1, info )
             ENDDO
             CALL pvmfsend( Xtid(-Hop(1)), msgtype, info ) 
           ENDIF

           IF ( direct .EQ. 2) THEN
             j  = Jlo(bo)
             i  = Ilo(bo)
             ir = Ihi(bo) - Ilo(bo) + 1
             DO k = k1, k2
               CALL pvmfpack( REAL8, A(i,j,k), ir, 1, info )
             ENDDO
             CALL pvmfsend( Ytid(-Hop(2)), msgtype, info ) 
           ENDIF

           IF ( direct .EQ. 3) THEN
             k  = Klo(bo)
             i  = Ilo(bo)
             ir = Ihi(bo) - Ilo(bo) + 1
             DO j = Jlo(bo), Jhi(bo)
               CALL pvmfpack( REAL8, A(i,j,k), ir, 1, info )
             ENDDO
             CALL pvmfsend( Ztid(-Hop(3)), msgtype, info ) 
           ENDIF

           ENDIF

         RETURN
       END

!----------------------------------------------------------------------
 
       SUBROUTINE Heard( A, bo, which )
!      --------------------------------
!PVM   Only used to receive interprocessor boundary talk
 
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         INTEGER bo, which, direct
         REAL A(0:N1,0:N2,0:N3)
         INTEGER Hop(3), i, j, k, k1, k2, ir, jr

!         IF (bo .GT. 6) RETURN ! only external boundaries

           k1 = Klo(bo)
           k2 = Khi(bo)

         DO k = 1,3
           Hop(k) = Orient(bo,k) ! step offsets
           IF (Hop(k) .NE. 0) direct = k
         ENDDO
 
!          Interprocessor talk b.c : receive
           IF ( Type(bo, which) .EQ. PVM_talk ) THEN
!            PRINT *, 'heard...', bo, which, Type(bo,which)
!            Receive A( i, j, k ) from next proc
             msgtype = which*10 + direct*Hop(direct) ! unique code
             
             IF ( direct .EQ. 1) THEN
               CALL pvmfrecv( Xtid(-Hop(1)), msgtype, info ) 
               j  = Jlo(bo)
               i  = Ilo(bo) -Hop(1)
               jr = Jhi(bo) - Jlo(bo) + 1
               DO k = k1, k2
                 CALL pvmfunpack( REAL8, A(i,j,k), jr, N2+1, info )
               ENDDO
             ENDIF

             IF ( direct .EQ. 2) THEN
               CALL pvmfrecv( Ytid(-Hop(2)), msgtype, info ) 
               j  = Jlo(bo) -Hop(2)
               i  = Ilo(bo)
               ir = Ihi(bo) - Ilo(bo) + 1
               DO k = k1, k2
                 CALL pvmfunpack( REAL8, A(i,j,k), ir, 1, info )
               ENDDO
             ENDIF
           
             IF ( direct .EQ. 3) THEN
               CALL pvmfrecv( Ztid(-Hop(3)), msgtype, info ) 
               k  = Klo(bo) -Hop(3)
               i  = Ilo(bo)
               ir = Ihi(bo) - Ilo(bo) + 1
               DO j = Jlo(bo), Jhi(bo)
                 CALL pvmfpack( REAL8, A(i,j,k), ir, 1, info )
               ENDDO
             ENDIF

           ENDIF

         RETURN
       END

!----------------------------------------------------------------------

       SUBROUTINE Press_Planes( P, Nx, Ny, Nz, k )
!      -------------------------------------------
!        Sets 2D guard planes on pressure, for plane k only
 
         INCLUDE 'jet3d.inc'
         INTEGER Nx, Ny, Nz, k, bound, oppos
         REAL P(0:Nx,0:Ny,0:Nz)

         DO bound = 1, 4 ! external only at present
           oppos = bound - (-1)**bound ! goes 2, 1, 4, 3
           CALL Guard_P( P, Nx, Ny, Nz, bound, 5, k)
           CALL Heard_P( P, Nx, Ny, Nz, oppos, 5, k)
         END DO

         RETURN
       END

!----------------------------------------------------------------------
 
       SUBROUTINE Guard_P( A, Nx, Ny, Nz, bo, which, k )
!      -------------------------------------------------
!      Sets 2D b.c. on pressure (in k3 space) for plane k
 
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         INTEGER Nx, Ny, Nz, bo, which
         REAL A(0:Nx,0:Ny,0:Nz)
         INTEGER i, j, k, i1, i2, j1, j2, ir, jr, direct

! messy --- can't we do better than this????         
         ir = Nx-1
         jr = Ny-1
         IF (bo .EQ. 1) THEN
           i1 = 1
           i2 = 0
           direct = 1
         ENDIF
         IF (bo .EQ. 2) THEN
           i1 = ir
           i2 = Nx
           direct = 1
         ENDIF
         IF (bo .EQ. 3) THEN
           j1 = 1
           j2 = 0
           direct = 2
         ENDIF
         IF (bo .EQ. 4) THEN
           j1 = jr
           j2 = Ny
           direct = 2
         ENDIF

         IF ( Type(bo, which) .EQ. Gradient ) THEN
           IF ( direct .EQ. 1 ) THEN
             DO j =  0, jr+1                     !orig. 1,jr
               A( i2, j, k ) = A( i1, j, k )
             ENDDO
           ELSE
             DO i =  0, ir+1                     !orig. 1,ir
               A( i, j2, k ) = A( i, j1, k )
             ENDDO
           ENDIF
         ENDIF
 
!          Interprocessor talk b.c : sending data
           IF ( Type(bo, which) .EQ. PVM_talk ) THEN
!            Send A( i, j, k ) to next proc
!            PRINT *, 'PVM called for ', bo, which
             CALL pvmfinitsend( PVMDEFAULT, info )

             IF ( direct .EQ. 1) THEN
               msgtype = which*10 + direct*(i1-i2) ! unique code
               CALL pvmfpack( REAL8, A(i1,1,k), jr, Ny+1, info )
               CALL pvmfsend( Xtid(i1-i2), msgtype, info ) 
             ELSE
               msgtype = which*10 + direct*(j1-j2) ! unique code
               CALL pvmfpack( REAL8, A(1,j1,k), ir, 1, info )
               CALL pvmfsend( Ytid(j1-j2), msgtype, info ) 
             ENDIF

           ENDIF

         RETURN
       END
 
!----------------------------------------------------------------------
 
       SUBROUTINE Heard_P( A, Nx, Ny, Nz, bo, which, k )
!      -------------------------------------------------
!PVM   Only used to receive interprocessor boundary talk
 
         INCLUDE 'jet3d.inc'
         INCLUDE 'pvm.inc'
         INTEGER Nx, Ny, Nz, bo, which, direct
         REAL A(0:Nx,0:Ny,0:Nz)
         INTEGER k, i1, i2, j1, j2, ir, jr
 
!         IF (bo .GT. 6) RETURN ! only external boundaries 

!        Interprocessor talk b.c : receive
         IF ( Type(bo, which) .EQ. PVM_talk ) THEN
!          Receive A( i, j, k ) from next proc

           ir = Nx-1         
           jr = Ny-1
! messy --- can't we do better than this????         
         ir = Nx-1
         jr = Ny-1
         IF (bo .EQ. 1) THEN
           i1 = 1
           i2 = 0
           direct = 1
         ENDIF
         IF (bo .EQ. 2) THEN
           i1 = ir
           i2 = Nx
           direct = 1
         ENDIF
         IF (bo .EQ. 3) THEN
           j1 = 1
           j2 = 0
           direct = 2
         ENDIF
         IF (bo .EQ. 4) THEN
           j1 = jr
           j2 = Ny
           direct = 2
         ENDIF

         IF ( direct .EQ. 1) THEN
           msgtype = which*10 + direct*(i2-i1) ! unique code
           CALL pvmfrecv( Xtid(i2-i1), msgtype, info ) 
           CALL pvmfunpack( REAL8, A(i2,1,k), jr, Ny+1, info )
         ENDIF

         IF ( direct .EQ. 2) THEN
           msgtype = which*10 + direct*(j2-j1) ! unique code
           CALL pvmfrecv( Ytid(j2-j1), msgtype, info ) 
           CALL pvmfunpack( REAL8, A(1,j2,k), ir, 1, info )
         ENDIF
           
         ENDIF

         RETURN
       END

!----------------------------------------------------------------------
      SUBROUTINE Read_inflow
!     ----------------------

      INCLUDE 'jet3d.inc'
      INTEGER i,j,k,idum,jdum

!     read boundary conditions from disk (as data file, ascii)   
      OPEN (UNIT=3, FILE='inflowbc.3d')
!      OPEN (UNIT=4, FILE='tempinflow.3d')
!     #changed
!     boundary 1
!     1. LEFT INFLOW PLANE X+ faced from 0 to N2
      DO j=0,N2
            READ (3,*) jdum, U1jet_in(j), U2jet_in(j)
      IF( j .NE. jdum )THEN
      write(6,*) ' Incompatible jet inflow data ! '
      STOP
      ENDIF
!
      ENDDO           
!
!     
!     2. TOP & BOTTOM INFLOW PLANE Y+- faced   
      DO i=1,M1
            READ (3,*) idum, U2en_top(i), U2en_bot(i)
      IF( i .NE. idum )THEN
      write(6,*) ' Incompatible entrainment inflow data !'
      STOP
      ENDIF
!
      ENDDO

!     3. TEMPERATURE INFLOW (LEFT PLANE ONLY)
!      DO j=0,M2
!            READ (4,*) jdum, T_in(j)
!      IF( j .NE. jdum )THEN
!      write(6,*) ' Incompatible jet inflow data ! '
!      STOP
!      ENDIF
!      ENDDO

!     
CLOSE(UNIT=4,STATUS='KEEP')

DO K=0, N3
   DO J=0, N2
      if( j .lt. jsplit) then
      T1(1,j,k) = 0.0
!      pass1(1,j,k) = T1(1,j,k)
      else
      T1(1,j,k) = 1.0
      endif
   ENDDO
ENDDO

DO I=1, M1
DO K=0, M3
T1(i,N2,k) = T1(i,M2,k)
T1(i,0,k) = T1(i,1,k)
ENDDO
ENDDO

!DO j=0,N2
!WRITE(*,*) j, T1(1,j,1)
!ENDDO

        
!     
!      CLOSE (UNIT=3,STATUS='KEEP')

!
      RETURN
      END
!-----------------------------------------------------------------------
 
       SUBROUTINE Inflow_bc(by,what)
!      ------------------------------
!      set jet inflow and entrainment b.c.

         INCLUDE 'jet3d.inc'
         INTEGER i,j,k,by,what,delay 
!        #changed
!        Boundary 1 or 7
!        1.LEFT INFLOW PLANE, left X+ faced
         IF (by .EQ. 1 .OR. by .EQ. 7) THEN     
!        U-FILE
         IF (what .EQ. 1) THEN
!        # jet inlet profile: ! bc 1 only
         DO j=0,N2                       ! 0 to M2+1
            DO k=0,N3
            IF( j.GE. 1 .AND. j.LE. M2 )THEN
            U1(1,j,k)=U1jet_in(j) + U1_ran(j,k)     
!            T1(i,j,k)= T_in(j)
            ELSE
            U1(1,j,k)=U1jet_in(j)
            ENDIF
            ENDDO
         ENDDO
         ENDIF
!        V-FILE
         IF (what .EQ. 2) THEN
         DO j=0,N2                       ! 0 to M2+1
            DO k=0,N3
            IF( j.GE. 1  .AND. j.LE. M2  )THEN
            U2(0,j,k)= U2jet_in(j) + U2_ran(j,k)   
            ELSE                      ! 
            U2(0,j,k)=U2jet_in(j)
            ENDIF
            ENDDO
         ENDDO

         ENDIF
!        W-FILE
         IF (what .EQ. 3) THEN 
         DO j=0,N2
            DO k=0,N3
            IF( j.GE. 1 .AND. j.LE. M2 )THEN
            U3(0,j,k)=U3_ran(j,k)
            ELSE
            U3(0,j,k)=0.0
            ENDIF
            ENDDO
         ENDDO
         ENDIF
         ENDIF      !by
!        2.TOP INFLOW PLANE, down Y- faced
         IF (by .EQ. 4) THEN
!        U-FILE
         IF (what .EQ. 1) THEN
         DO i=1,M1
            DO k=0,N3
                  U1(i,N2,k)=0.0
            ENDDO
         ENDDO
         ENDIF
!        V-FILE
         IF (what .EQ. 2) THEN
         DO i=1,M1
            DO k=0,N3
                  U2(i,N2,k)=U2en_top(i)
            ENDDO
         ENDDO
         ENDIF
!        W-FILE
         IF (what .EQ. 3) THEN
         DO i=1,M1
            DO k=0,N3
                  U3(i,N2,k)=0.0
            ENDDO
         ENDDO
         ENDIF
      ENDIF !by

!
!        3.BOTTOM INFLOW PLANE, up Y+ faced
         IF (by .EQ. 3) THEN
!        U-FILE  j=0
         IF (what .EQ. 1) THEN
         DO i=1,M1
            DO k=0,N3
                     U1(i,0,k)=0.0
            ENDDO
         ENDDO
         ENDIF
!        V-FILE    j=1
         IF (what .EQ. 2) THEN
         DO i=1,M1
            DO k=0,N3
                     U2(i,1,k)=U2en_bot(i)
            ENDDO
         ENDDO
         ENDIF
!
!        W-FILE    j=0
         IF (what .EQ.  3) THEN
         DO i=1,M1
            DO k=0,N3
                     U3(i,0,k)=0.0
            ENDDO
         ENDDO
         ENDIF

         ENDIF    !by
!
!

         IF (by .EQ. 5) THEN
!        T-FILE  k=0
         IF (what .EQ. 4) THEN
            DO i=1,M1
            DO j=0,N2
                     T1(i,j,0)= T1(i,j,M3)
            ENDDO
            ENDDO
            ENDIF
            ENDIF

         IF (by .EQ. 6) THEN
!        T-FILE  k=0
         IF (what .EQ. 4) THEN
            DO i=1,M1
            DO j=0,N2
                     T1(i,j,N3)= T1(i,j,1)
            ENDDO
            ENDDO
            ENDIF
            ENDIF

         RETURN
       END
 
!----------------------------------------------------------------------
      SUBROUTINE Read_outflow
!     -----------------------
!     read outflow condition to define convective velocity    
!     # NOT USED.
      INCLUDE 'jet3d.inc'
      INTEGER j
!      OPEN (UNIT=3, FILE='outflowc.3d', FORM='UNFORMATTED') 
!            READ (3) SUMOUTBL
!            READ (3) (U1bl(j),j=0,N2)
!      CLOSE (UNIT=3, STATUS='KEEP')

      RETURN
      END
!---------------------------------------------------------------------
 
       SUBROUTINE Outflow_bc
!      ---------------------
!      set outflow condition

         INCLUDE 'jet3d.inc'
         INTEGER I, J, K, bs
         REAL SUMIN, SUMOUT, UMout, U1_out (0:N2), sum

         SUMIN = 0.0
!        ******************** OLD VERSION OF COMP. OF INFLOW FLUX
!         DO J=Jlo(1),Jhi(1)
!            DO K=1,M3
!                  SUMIN=SUMIN+U1(1,J,K)*dy(J)
!            ENDDO
!         ENDDO
!        ******************** END OF OLD VERSION
!                
!        ******************** computation of inflow flux

!        Inflow Boundary faced in X-Direction
                  DO J = Jlo(1),Jhi(1)           !boundary 1 only
                        DO K = Klo(1), Khi(1)
                        SUMIN = SUMIN + U1(1,J,K)*dy(J)  
                        ENDDO
                  ENDDO
!          write(6,*) 'jet inflow flux =',SUMIN
!  
         
!        Inflow Boundary faced in Y-Direction, TOP & BOTTOM bound..
!        Upper boundary 4
                  DO I = 1, m1
                  DO K = Klo(4), Khi(4)
                  SUMIN = SUMIN + abs( U2(I,m2+1,K)*dx(I) )
                  ENDDO
                  ENDDO
!
!         Lower boundary 3
!         
                  DO I = 1, m1
                  DO K = Klo(3), Khi(3)
                  SUMIN = SUMIN + U2(I,1,K)*dx(I)
                  ENDDO
                  ENDDO
!
!
!        ******************** end of computation of inflow flux
!
!        ******************** outflow boundary condition 
!        2. setting convective outflow boundary conditions
!        OUTFLOW boundary is the right X-Plane (X- facing)

         UMout = SUMIN / ((y(Jhi(2)+1)-y(Jlo(2)))*M3)

!         DO J = 0, N2
! NOT USED           U1_out(J)=SUMIN*U1bl(J)/SUMOUTBL
!         ENDDO                                                    !to here

         IF (Step_no .EQ. 1 .AND. New_Run) THEN
          DO K = 0, N3
          DO J = 0, N2
            U1(N1,J,K) = UMout
            U2(N1,J,K) = 0
            U3(N1,J,K) = 0
          ENDDO
          ENDDO
         ELSE
!         Convective b.c.     CONSTANT  convection velocity !swap U1_out<>UMout
          DO K = 0, N3
          DO J = 0, N2
            U1(N1,J,K) = U1(N1,J,K) - Delta_t*UMout*               &
            (U1(N1,J,K)-U1(M1,J,K))/dx(M1)
            U2(N1,J,K) = U2(N1,J,K) - Delta_t*UMout*               &
            (U2(N1,J,K)-U2(M1,J,K))/(0.5*(dx(N1)+dx(M1)))
            U3(N1,J,K) = U3(N1,J,K) - Delta_t*UMout*               &
            (U3(N1,J,K)-U3(M1,J,K))/(0.5*(dx(M1)+dx(N1)))
          ENDDO
          ENDDO
         ENDIF
!        end of test ############################################
! 
          SUMOUT = 0.0    
          DO J = Jlo(2), Jhi(2)
          DO K = 1, M3
            SUMOUT = SUMOUT + U1(N1,J,K)*dy(J)
          ENDDO
          ENDDO
!         scaling velocities
          IF(ABS(SUMOUT) .GT. 0.0) THEN
            DO K = 0, N3
            DO J = 0, N2
              U1(N1,J,K) = SUMIN*U1(N1,J,K)/SUMOUT
              U2(N1,J,K) = SUMIN*U2(N1,J,K)/SUMOUT 
              U3(N1,J,K) = SUMIN*U3(N1,J,K)/SUMOUT 
            ENDDO
            ENDDO
          ENDIF
!         
!        ******************** end of outflow boundary condition
!         
!        ******************** computation of outflow flux
          SUMOUT = 0.0
          DO J = Jlo(2), Jhi(2)
          DO K = 1, M3
            SUMOUT = SUMOUT + U1(N1,J,K)*dy(J)
          ENDDO
          ENDDO
!        ******************** end of outflow flux


           WRITE (*,*) 'INFLOW/OUTFLOW Check at Time Step= ',Step_no
           WRITE (*,100) 'SUMIN= ',SUMIN, '  SUMOUT= ', SUMOUT,         &
                        '  UMout= ',UMout


100       FORMAT (a7,e17.10,a10,e17.10,a9,e17.10)
          RETURN
        END
!----------------------------------------------------------------------

      SUBROUTINE Outflow_temp
!      ---------------------
!      set outflow condition

         INCLUDE 'jet3d.inc'
         INTEGER I, J, K, bs
         REAL SUMIN, SUMOUT, UMout, U1_out (0:N2), sum

IF (Step_no .EQ. 1 .AND. New_Run) THEN
          DO K = 0, N3
          DO J = 0, N2
            T1(N1,J,K) = T1(1,j,k) 
          ENDDO
          ENDDO
         ELSE
!         Convective b.c.     CONSTANT  convection velocity !swap U1_out<>UMout
          DO K = 0, N3
          DO J = 0, N2
            T1(N1,J,K) = T1(M1,j,k)!
          ENDDO
          ENDDO
         ENDIF
!        end of test ############################################


RETURN
END SUBROUTINE Outflow_temp
!--------------------------------------------------------------------


        SUBROUTINE Ini_Ran
!       ------------------
!       set new random numbers 
        INCLUDE 'jet3d.inc'
        INTEGER k,j
!        REAL(KIND=2) :: RANDOM, RANVEC(N2,N3)
        REAL :: RANF
!
        DO k = 0, N3
        DO j = jsplit, N2

           CALL RANDOM_NUMBER(RANF)
        U1_ran(j,k)=0.005*(u1jet_in(m2-1))*2.0*(RANF - 0.5)
         CALL RANDOM_NUMBER(RANF)
        U2_ran(j,k)=0.005*(u1jet_in(m2-1))*2.0*(RANF- 0.5) !0.5%disturbances
         CALL RANDOM_NUMBER(RANF)
        U3_ran(j,k)=0.005*(u1jet_in(m2-1))*2.0*(RANF - 0.5)
        ENDDO
        ENDDO

        DO k = 0, N3
        DO j = 0, jsplit-1
            CALL RANDOM_NUMBER(RANF)
        U1_ran(j,k)=0.005*(u1jet_in(2))*2.0*(RANF - 0.5)
         CALL RANDOM_NUMBER(RANF)
        U2_ran(j,k)=0.005*(u1jet_in(2))*2.0*(RANF- 0.5) !0.5%disturbances
         CALL RANDOM_NUMBER(RANF)
        U3_ran(j,k)=0.005*(u1jet_in(2))*2.0*(RANF - 0.5)
        ENDDO
        ENDDO



        RETURN
        END
!

!----------------------------------------------------------------------

        SUBROUTINE Div_U (more,sum)
!       ---------------------------
!       Calculation of DIVERGENCE (Ui)   

        INCLUDE 'jet3d.inc'

        REAL sum,max,div
        INTEGER i,j,k,imax,jmax,kmax
        LOGICAL more

        sum=0.0
        max=0.0
        DO i = 1, M1
        DO j = 1, M2
        DO k = 1, M3
            div=(U1(i+1,j,k)-U1(i,j,k))/dx(i)+(U2(i,j+1,k)             &
            -U2(i,j,k))/dy(j)+(U3(i,j,k+1)-U3(i,j,k))/dz(k)
        sum=sum + div*dx(i)*dy(j)*dz(k)   !=INT div U dV
            IF (ABS(div) .GT. ABS(max)) THEN
                  max=div
                  imax=i
                  jmax=j
                  kmax=k
            ENDIF
        ENDDO
        ENDDO
        ENDDO
        IF (Debug .GT. 0) THEN
            WRITE (*,*) 'Continuity check    : '
            WRITE (*,*) 'INT div(u) dV ='
            WRITE (*,100) sum
            WRITE (*,*) 'Max. residual =',max,' at ',imax,jmax,kmax
        ENDIF
        IF (max .GE. accuracy) THEN
          more=.TRUE.
        ELSE 
          more=.FALSE.
        ENDIF
 
100     FORMAT (e17.10)
        RETURN
        END
!       ################### end of div_u ##################

       SUBROUTINE Monitor
!      ------------------
!
         INCLUDE 'jet3d.inc'
         REAL Co, CoX, CoY, CoZ, Tn, Xn, Yn, Zn, divergence
         INTEGER i, j, k, ip, jp, kp
         LOGICAL dummy
!
         CoX = 0.0
         CoY = 0.0
         CoZ = 0.0
         Co = 0.0
         DO k = 1, M3
         DO j = 1, M2
         DO i = 1, M1
           Xn = Delta_t * ABS( U1(i,j,k) / dx(i) )
           Yn = Delta_t * ABS( U2(i,j,k) / dy(j) )
           Zn = Delta_t * ABS( U3(i,j,k) / dz(k) )
           Tn = Xn + Yn + Zn
           IF (Tn .GT. Co ) THEN
             Cox = Xn
             Coy = Yn
             Coz = Zn
             Co  = Tn
             ip = i
             jp = j
             kp = k
           ENDIF
         ENDDO
         ENDDO
         ENDDO
!                          

       CPUtot_d=CPUtot_n-CPUtot_o
       CPUtot_o=CPUtot_n
       write(*,'(''*Info: CFL number: '',g15.8, '' for time step '',i10)') Co,Step_No

!       WRITE (9,198) Step_No, Time, Co, Cox, Coy, Coz, i, j, k, divergence, CPUtot_d, CPUpr_d

!
!198    FORMAT (i6,e14.6,e14.6,e14.6)
!198    FORMAT (i6,e14.6,e14.6,e14.6,e14.6,e14.6,i6,i6,i6,e14.6,e14.6,e14.6)
!
       RETURN
       END

!----------------------------------------------------------------------

      SUBROUTINE Field_out
!     --------------------
        INCLUDE 'jet3d.inc'
        INTEGER i, j, k
        CHARACTER*1 T_index(6)
        CHARACTER*12 field
        REAL Ev, sum
        REAL :: xr(1:M1,1:M2), yr(1:M1,1:M2), zr(1:M1,1:M2), tr(1:M1,1:M2)
        REAL :: xrz(1:M1,1:M3), yrz(1:M2,1:M3), zrz(1:M1,1:M3), trz(1:M1,1:M3)
        REAL :: tempt(1:M1,1:M2), tempu(1:M1,1:M2), tempv(1:M1,1:M2), tempw(1:M1,1:M2)
        real :: tempa(1:M1,1:M2)

       CALL Int2char(Step_No,T_index)
       field='r2'//T_index(1)//T_index(2)//T_index(3)//T_index(4)      &
       //T_index(5)//T_index(6)//'.3d'



       WRITE(6,*) 'Writing out flow visualsiation files'
       WRITE(99,*) 'VARIABLES = "x (m)" "y (m)" "u1" "u2" "u3" "ev"'
       WRITE(99,*) 'ZONE I=', m2,', J=', m1
       WRITE(90,*) 'VARIABLES = "x (m)" "y (m)" "T" "pass1" "NO" "O3" "NO2"'
       WRITE(90,*) 'ZONE I=', m2,', J=', m1
       DO i=1,m1
       DO j=1,m2
       if( m3 .eq. 1) then
       WRITE(90,*) x(i),y(j), T1(i,j,1), pass1(i,j,1), rNO(i,j,1),rO3(i,j,1), rNO2(i,j,1)
       Ev= ( E(i,j,1)-Nu )*100/Nu
       WRITE(99,*) x(i),y(j),U1(i,j,1),U2(i,j,1),U3(i,j,1),Ev
       else
       WRITE(90,*) x(i),y(j), T1(i,j,m3/2), pass1(i,j,m3/2), rNO(i,j,m3/2),rO3(i,j,m3/2), rNO2(i,j,m3/2)
       Ev= ( E(i,j,m3/2)-Nu )*100/Nu
       WRITE(99,*) x(i),y(j),U1(i,j,m3/2),U2(i,j,m3/2),U3(i,j,m3/2),Ev
       endif
       ENDDO
       ENDDO
       
!work out spanwise averaged quantities

         DO i=1,M1
            DO j=1,M2
                  sum=0.0
                  DO k=1,M3
                  sum=sum+T1(i,j,k)*dz(k)
                  ENDDO !k
            tempt(i,j)=sum/(z(M3+1)-z(1))
            ENDDO
      ENDDO   
      
               DO i=1,M1
            DO j=1,M2
                  sum=0.0
                  DO k=1,M3
                  sum=sum+pass1(i,j,k)*dz(k)
                  ENDDO !k
            tempu(i,j)=sum/(z(M3+1)-z(1))
            ENDDO
      ENDDO
      
DO i=1,M1
DO j=1,M2
sum=0.0
DO k=1,M3
sum=sum+rNO(i,j,k)*dz(k)
ENDDO !k
tempv(i,j)=sum/(z(M3+1)-z(1))
ENDDO
ENDDO

DO i=1,M1
DO j=1,M2
sum=0.0
DO k=1,M3
sum=sum+rO3(i,j,k)*dz(k)
ENDDO !k
tempw(i,j)=sum/(z(M3+1)-z(1))
ENDDO
ENDDO

DO i=1,M1
DO j=1,M2
sum=0.0
DO k=1,M3
sum=sum+rNO2(i,j,k)*dz(k)
ENDDO !k
tempa(i,j)=sum/(z(M3+1)-z(1))
ENDDO
ENDDO

      WRITE(91,*) 'VARIABLES = "x (m)" "y (m)" "T ave" "pass1" "NO" "O3" "NO2"'
       WRITE(91,*) 'ZONE I=', M2,', J=', M1
        DO I=1, M1
        DO J=1,M2
        WRITE(91,*) x(i), y(j), tempt(i,j), tempu(i,j), tempv(i,j), tempw(i,j), tempa(i,j)
        ENDDO
        ENDDO

!Calculate the vorticity field in the flow

      !x-component
      DO I=1, M1
         DO J=1, M2
            DO K=1, M3
               Wa(i,j,k) = -(U2(i,j,k) - U2(i,j,k-1))/(z(k) - z(k-1)) + &
            ((U3(i,j,k) - U3(i,j-1,k)) / (y(j) - y(j-1)))
            ENDDO
         ENDDO
      ENDDO
     
      
      DO i=1,M1
            DO j=1,M2
                  sum=0.0
                  DO k=1,M3
                  sum=sum+Wa(i,j,k)*dz(k)
                  ENDDO !k
            xr(i,j)=sum/(z(M3+1)-z(1))
            ENDDO
      ENDDO

      !y-component
      DO I=1, M1
         DO J=1, M2
            DO K=1, M3
               Wb(i,j,k) = -(U3(i,j,k) - U3(i-1,j,k))/(x(i) - x(i-1)) + &
            ((U1(i,j,k) - U1(i,j,k-1)) / (z(k) - z(k-1)))
            ENDDO
         ENDDO
      ENDDO
     
      
       DO i=1,M1
            DO j=1,M2
                  sum=0.0
                  DO k=1,M3
                  sum=sum+Wb(i,j,k)*dz(k)
                  ENDDO !k
            yr(i,j)=sum/(z(M3+1)-z(1))
            ENDDO
      ENDDO

      ! z-component
      DO I=1, M1
         DO J=1, M2
            DO K=1, M3
               Wc(i,j,k) = (U2(i,j,k) - U2(i-1,j,k))/(x(i) - x(i-1)) - &
            ((U1(i,j,k) - U1(i,j-1,k)) / (y(j) - y(j-1)))
            ENDDO
         ENDDO
      ENDDO

      
    
      
      DO i=1,M1
            DO j=1,M2
                  sum=0.0
                  DO k=1,M3
                  sum=sum+Wc(i,j,k)*dz(k)
                  ENDDO !k
            zr(i,j)=sum/(z(M3+1)-z(1))
            ENDDO
      ENDDO

    DO I=1,M1
    DO J=1, M2
    tr(i,j) = xr(i,j) + yr(i,j) + zr(i,j)
    ENDDO
    ENDDO

      !Write outputs to file
      WRITE(65,*) 'VARIABLES = "x" "y" "x-vort" "z-vort" "t-vort"'
      WRITE(65,*) 'ZONE I= ', M2, ', J=', M1
      DO I=1, M1
         DO J=1, M2
            WRITE(65,*) x(i), y(j), xr(i,j), zr(i,j), tr(i,j)
         ENDDO
      ENDDO
      
!output to a y-z plane, pre and post transition
      WRITE(94,*) 'VARIABLES = "x" "y" "v-vel" "u-vel" "z-vort" "w-vel" "temp"'
      WRITE(94,*) 'ZONE I= ', M3, ', J=', M2
      DO J=1, M2
         DO K=1, M3
            WRITE(94,*) y(j), z(k), U2(166,j,k), U1(166,j,k), Wc(166,j,k), U2(166,j,k), T1(166,j,k)
         ENDDO
      ENDDO
      
 

197    FORMAT (a11,i6,a16,f12.6) 
200    FORMAT (e19.10,e19.10,e19.10,e19.10)
201    FORMAT (a11)

        RETURN
      END

!----------------------------------------------------------------------

       SUBROUTINE Dump
!      ---------------

         INCLUDE 'jet3d.inc'
         INTEGER i, j, k
         character*3 int_as_String
         character*12 all

         write(int_as_string, fmt='(I3.3)') myrank
         all='sfield'//int_as_string//'.3d'

!
!        1. write old velocity, pressure & history fields 
!        performed only on master
         if(myrank .eq. 0) then
         OPEN (UNIT=8,FILE='dump.3d',FORM='UNFORMATTED')
         REWIND(8)
         WRITE (8) Step_No, Time
         WRITE (8) N1,N2,N3
         WRITE (8) LES
         WRITE (8) (((U1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((U2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((U3(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((P1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((H1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((H2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((H3(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((T1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         WRITE (8) (((pass1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write (8) (((rNO(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write (8) (((rO3(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write (8) (((rNO2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write (8) (((rO2(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         CLOSE (UNIT=8, STATUS='KEEP')
         endif

!        2.  Write stochastic fields from all cores
         open(unit=1000+myrank,form='unformatted')
         rewind(1000+myrank)
         write(1000+myrank) (((pasf1(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write(1000+myrank) (((NOf(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write(1000+myrank) (((O3f(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write(1000+myrank) (((NO2f(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         write(1000+myrank) (((O2f(i,j,k),i=0,N1),j=0,N2),k=0,N3)
         close(unit=1000+myrank, status='keep')

        !Changed this to make output an ASCII file
if(myrank .eq. 0) then
OPEN ( UNIT=8, FILE='stat.3d')
REWIND(8)
WRITE(8,*) No_Sam
DO I=1,M1
DO J=1,M2
WRITE(8,*) u1_acc(i,j),u2_acc(i,j),u3_acc(i,j),p1_acc(i,j),T1_acc(i,j)
ENDDO
ENDDO

DO I=1,M1
DO J=1,M2
WRITE(8,*) uu_acc(i,j),vv_acc(i,j),ww_acc(i,j),pp_acc(i,j),vp_acc(i,j)
ENDDO
ENDDO

DO I=1,M1
DO J=1,M2
WRITE(8,*) uv_acc(i,j),uw_acc(i,j),vw_acc(i,j),up_acc(i,j),wp_acc(i,j)
ENDDO
ENDDO

DO I=1,M1
DO J=1,M2
WRITE(8,*) cx_acc(i,j),cy_acc(i,j),cz_acc(i,j),tt_acc(i,j),T1_flu(i,j)
ENDDO
ENDDO

DO I=1,M1
DO J=1,M2
WRITE(8,*) uu_flu(i,j),vv_flu(i,j),ww_flu(i,j),uv_flu(i,j),uw_flu(i,j),vw_flu(i,j)
ENDDO
ENDDO

DO I=1,M1
DO J=1,M2
WRITE(8,*) pass1_acc(i,j), pp1_acc(i,j),pass1_flu(i,j)
ENDDO
ENDDO

DO I=1, M1
DO J=1, M2
WRITE(8,*) rNO_acc(i,j), ppNO_acc(i,j),rNO_flu(i,j)
ENDDO
ENDDO

DO I=1, M1
DO J=1, M2
WRITE(8,*) rO3_acc(i,j), ppO3_acc(i,j),rO3_flu(i,j)
ENDDO
ENDDO

DO I=1, M1
DO J=1, M2
WRITE(8,*) rNO2_acc(i,j), ppNO2_acc(i,j),rNO2_flu(i,j)
ENDDO
ENDDO

DO I=1, M1
DO J=1, M2
WRITE(8,*) rO2_acc(i,j), ppO2_acc(i,j),rO2_flu(i,j)
ENDDO
ENDDO
endif !myrank


CLOSE (UNIT=8,STATUS='KEEP')

         RETURN

       END

!-------------------------------------------------------------------

       SUBROUTINE Statistics
!     ---------------------
      INCLUDE 'jet3d.inc'
!
      INTEGER i,j,k
      REAL u_cc,v_cc,w_cc,curlx_cc,curly_cc,curlz_cc
      REAL red1(1:M1,1:M2),red2(1:M1,1:M2),red3(1:M1,1:M2),            &
           red4(1:M1,1:M2),red5(1:M1,1:M2)
!     temporary work fields to save interpolated values
!
!     U, P - STATISTICS
!     1. interpolate velocities to the cell center, U1 saved in Wa, U2
!     in Wb & U3 in Wc
      DO i=1,M1
      DO j=1,M2
      DO k=1,M3
      Wa(i,j,k)=0.5*(U1(i,j,k)+U1(i+1,j,k))
      Wb(i,j,k)=0.5*(U2(i,j,k)+U2(i,j+1,k))
      Wc(i,j,k)=0.5*(U3(i,j,k)+U3(i,j,k+1))
      ENDDO
      ENDDO
      ENDDO
!
!     2. reduction of u1, u2, u3, p1 in spanwise direction 
!     and accumulation
      CALL red_z (Wa,red1,1,M1,1,M2,1,M3)  !U1-Average
      CALL red_z (Wb,red2,1,M1,1,M2,1,M3)  !U2-Average
      CALL red_z (Wc,red3,1,M1,1,M2,1,M3)  !U3-Average
      CALL red_z (P1,red4,1,M1,1,M2,1,M3)  !P1-Average 
      CALL red_z (T1,red5,1,M1,1,M2,1,M3)  !T1-Average
      DO i=1,M1
      DO j=1,M2
      u1_acc(i,j)=(u1_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      u2_acc(i,j)=(u2_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      u3_acc(i,j)=(u3_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      p1_acc(i,j)=(p1_acc(i,j)*(No_Sam-1)+red4(i,j))/No_Sam
      T1_acc(i,j)=(T1_acc(i,j)*(No_Sam-1)+red5(i,j))/No_Sam
      ENDDO
      ENDDO

      !fdf statistics
      CALL red_z (pass1,red1,1,M1,1,M2,1,M3)  !pass1-Average
      CALL red_z (rNO,red2,1,M1,1,M2,1,M3)
      CALL red_z (rO3,red3,1,M1,1,M2,1,M3)
      CALL red_z (rNO2,red4,1,M1,1,M2,1,M3)
      CALL red_z (rO2,red5,1,M1,1,M2,1,M3)
      DO I=1, M1
      DO J=1, M2 
      pass1_acc(i,j)=(pass1_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      rNO_acc(i,j)=(rNO_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      rO3_acc(i,j)=(rO3_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      rNO2_acc(i,j)=(rNO2_acc(i,j)*(No_Sam-1)+red4(i,j))/No_Sam
      rO2_acc(i,j)=(rO2_acc(i,j)*(No_Sam-1)+red5(i,j))/No_Sam
      ENDDO
      ENDDO

      
!
!     UU, VV, WW, PP STATISTICS
!     1. computation of moments, uu saved in Wa, vv in Wb, ww in Wc
!     uses center cell velocities from above, pp saved in Wd
      DO i=1,M1
      DO j=1,M2
      DO k=1,M3
      Wa(i,j,k)=Wa(i,j,k)*Wa(i,j,k)
      Wb(i,j,k)=Wb(i,j,k)*Wb(i,j,k)
      Wc(i,j,k)=Wc(i,j,k)*Wc(i,j,k)
      Wd(i,j,k)=P1(i,j,k)*P1(i,j,k) 
      We(i,j,k)=T1(i,j,k)*T1(i,j,k)
      Wg(i,j,k)=pass1(i,j,k)*pass1(i,j,k)
      ENDDO
      ENDDO
      ENDDO
!
!     2. reduction of uu, vv, ww, pp in spanwise direction 
!     and accumulation
      CALL red_z (Wa,red1,1,M1,1,M2,1,M3)  !U1U1-Average
      CALL red_z (Wb,red2,1,M1,1,M2,1,M3)  !U2U2-Average
      CALL red_z (Wc,red3,1,M1,1,M2,1,M3)  !U3U3-Average
      CALL red_z (Wd,red4,1,M1,1,M2,1,M3)  !P1P1-Average
      CALL red_z (We,red5,1,M1,1,M2,1,M3)  !T1T1-Average
      DO i=1,M1
      DO j=1,M2
      uu_acc(i,j)=(uu_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      vv_acc(i,j)=(vv_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      ww_acc(i,j)=(ww_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      pp_acc(i,j)=(pp_acc(i,j)*(No_Sam-1)+red4(i,j))/No_Sam  
      tt_acc(i,j)=(tt_acc(i,j)*(No_Sam-1)+red5(i,j))/No_Sam
      ENDDO
      ENDDO



      DO i=1,M1
      DO j=1,M2
      DO k=1,M3
      Wa(i,j,k)=pass1(i,j,k)*pass1(i,j,k)
      Wb(i,j,k)=rNO(i,j,k)*rNO(i,j,k)
      Wc(i,j,k)=rO3(i,j,k)*rO3(i,j,k)
      Wd(i,j,k)=rNO2(i,j,k)*rNO2(i,j,k)
      We(i,j,k)=rO2(i,j,k)*rO2(i,j,k)
      ENDDO
      ENDDO
      ENDDO

      CALL red_z (Wa,red1,1,M1,1,M2,1,M3)  !pass1**2-Average
      CALL red_z (Wb,red2,1,M1,1,M2,1,M3)
      CALL red_z (Wc,red3,1,M1,1,M2,1,M3)
      CALL red_z (Wd,red4,1,M1,1,M2,1,M3)
      CALL red_z (We,red5,1,M1,1,M2,1,M3)
      DO I=1, M1
      DO J=1, M2
      pp1_acc(i,j)=(pp1_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      ppNO_acc(i,j)=(ppNO_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      ppO3_acc(i,j)=(ppO3_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      ppNO2_acc(i,j)=(ppNO2_acc(i,j)*(No_Sam-1)+red4(i,j))/No_Sam
      ppO2_acc(i,j)=(ppO2_acc(i,j)*(No_Sam-1)+red5(i,j))/No_Sam
      ENDDO
      ENDDO
      

      !fluctuation calcs
      DO I = 1, M1
         DO J = 1, M2
            uu_flu(I,J) = ( uu_flu(i,j)*(No_Sam-1) + (uu_acc(I,J) - u1_acc(I,J)**2) )/No_Sam
            vv_flu(I,J) = ( vv_flu(i,j)*(No_Sam-1) + (vv_acc(I,J) - u2_acc(I,J)**2) )/No_Sam
            ww_flu(I,J) = ( ww_flu(i,j)*(No_Sam-1) + (ww_acc(I,J) - u3_acc(I,J)**2) )/No_Sam
            T1_flu(i,j) = ( T1_flu(i,j)*(No_Sam-1) + (tt_acc(I,J) - t1_acc(I,J)**2) )/No_Sam
            pass1_flu(i,j) = ( pass1_flu(i,j)*(No_Sam-1) + (pp1_acc(I,J) - pass1_acc(I,J)**2) )/No_Sam
            rNO_flu(I,J) = ( rNO_flu(i,j)*(No_Sam-1) + (ppNO_acc(I,J) - rNO_acc(I,J)**2) )/No_Sam
            rO3_flu(I,J) = ( rO3_flu(i,j)*(No_Sam-1) + (ppO3_acc(I,J) - rO3_acc(I,J)**2) )/No_Sam
            rNO2_flu(I,J) = ( rNO2_flu(i,j)*(No_Sam-1) + (ppNO2_acc(I,J) - rNO2_acc(I,J)**2) )/No_Sam
            rO2_flu(I,J) = ( rO2_flu(i,j)*(No_Sam-1) + (ppO2_acc(I,J) - rO2_acc(I,J)**2) )/No_Sam
          ENDDO
      ENDDO

!
!     UV, UW, VW STATISTICS
!     1. computation of moments, uv saved in Wa, uw in Wb, vw in Wc
      DO i=1,M1
      DO j=1,M2
      DO k=1,M3
      Wa(i,j,k)=0.25*(U1(i,j,k)+U1(i+1,j,k))*(U2(i,j,k)+U2(i,j+1,k)) 
      Wb(i,j,k)=0.25*(U1(i,j,k)+U1(i+1,j,k))*(U3(i,j,k)+U3(i,j,k+1))
      Wc(i,j,k)=0.25*(U2(i,j,k)+U2(i,j+1,k))*(U3(i,j,k)+U3(i,j,k+1))
      ENDDO
      ENDDO
      ENDDO
!
!     2. reduction of uv, uw, vw in spanwise direction 
!     and accumulation
      CALL red_z (Wa,red1,1,M1,1,M2,1,M3)  !U1U2-Average
      CALL red_z (Wb,red2,1,M1,1,M2,1,M3)  !U1U3-Average
      CALL red_z (Wc,red3,1,M1,1,M2,1,M3)  !U2U3-Average
      DO i=1,M1
      DO j=1,M2
      uv_acc(i,j)=(uv_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      uw_acc(i,j)=(uw_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      vw_acc(i,j)=(vw_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      ENDDO
      ENDDO
      
      !Fluctuations
      DO I = 1, M1
         DO J = 1, M2
            uv_flu(I,J) = ( uv_flu(i,j)*(No_Sam-1) + (uv_acc(I,J) - u1_acc(I,J)**2) )/No_Sam
            uw_flu(I,J) = ( uw_flu(i,j)*(No_Sam-1) + (uw_acc(I,J) - u1_acc(I,J)**2) )/No_Sam
            vw_flu(I,J) = ( vw_flu(i,j)*(No_Sam-1) + (vw_acc(I,J) - u2_acc(I,J)**2) )/No_Sam
         ENDDO
      ENDDO

!     UP, VP, WP STATISTICS
!     1. computation of moments, up saved in Wa, vp in Wb, wp in Wc;
      DO i=1,M1
      DO j=1,M2
      DO k=1,M3
      Wa(i,j,k)=0.5*(U1(i,j,k)+U1(i+1,j,k))*P1(i,j,k) 
      Wb(i,j,k)=0.5*(U2(i,j,k)+U2(i,j+1,k))*P1(i,j,k)
      Wc(i,j,k)=0.5*(U3(i,j,k)+U3(i,j,k+1))*P1(i,j,k)
      ENDDO
      ENDDO
      ENDDO
!
!     2. reduction of up, vp, wp in spanwise direction 
!     and accumulation
      CALL red_z (Wa,red1,1,M1,1,M2,1,M3)  !U1P1-Average
      CALL red_z (Wb,red2,1,M1,1,M2,1,M3)  !U2P1-Average
      CALL red_z (Wc,red3,1,M1,1,M2,1,M3)  !U3P1-Average
      DO i=1,M1
      DO j=1,M2
      up_acc(i,j)=(up_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      vp_acc(i,j)=(vp_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      wp_acc(i,j)=(wp_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      ENDDO
      ENDDO
!
!     U x CURL U  STATISTICS
!     1. compute rot_u 
!     compute curl(u) on an orthogonal 3D staggered grid
!     (curl)z-comp. is defined at SW corner of z-plane trough the cell c.
!     (curl)y-comp. is ....... at SW corner of y-plane ....
!     (curl)x-comp. is ....... at SW corner of x-plane ....
!
!     x-component on work arr. Wa
!     y-component on work arr. Wb 
!     z-component on work arr. Wc 
      DO i=1,N1
      DO j=1,N2
      DO k=1,N3
            Wa(i,j,k)=(U3(i,j,k)-U3(i,j-1,k))*dy1(j)                   &
                     -(U2(i,j,k)-U2(i,j,k-1))*Rdz
            Wb(i,j,k)=(U1(i,j,k)-U1(i,j,k-1))*Rdz                      &
                     -(U3(i,j,k)-U3(i-1,j,k))*dx1(i)
            Wc(i,j,k)=(U2(i,j,k)-U2(i-1,j,k))*dx1(i)                   &
                     -(U1(i,j,k)-U1(i,j-1,k))*dy1(j)
      ENDDO
      ENDDO
      ENDDO
!
!     2. computation of u x curl_u components at the cell center
!     x-component saved in Wa, y-component saved in Wb, z-component
!     saved in Wc
!     NOTE: not recursive as forward looking!
      DO i=1,M1
      DO j=1,M2
      DO k=1,M3
      u_cc=0.5*(U1(i,j,k)+U1(i+1,j,k))
      v_cc=0.5*(U2(i,j,k)+U2(i,j+1,k))
      w_cc=0.5*(U3(i,j,k)+U3(i,j,k+1))
      curlx_cc=0.25*(Wa(i,j,k)+Wa(i,j+1,k)+Wa(i,j+1,k+1)+Wa(i,j,k+1))
      curly_cc=0.25*(Wb(i,j,k)+Wb(i,j,k+1)+Wb(i+1,j,k+1)+Wb(i+1,j,k))
      curlz_cc=0.25*(Wc(i,j,k)+Wc(i,j+1,k)+Wc(i+1,j+1,k)+Wc(i+1,j,k))
!
      Wa(i,j,k)=v_cc*curlz_cc - w_cc*curly_cc   !x-comp. of U x rotU
      Wb(i,j,k)=w_cc*curlx_cc - u_cc*curlz_cc   !y-comp. of U x rotU
      Wc(i,j,k)=u_cc*curly_cc - v_cc*curlx_cc   !z-comp. of U x rotU 
      ENDDO
      ENDDO
      ENDDO
!
!     3. reduction of U x rotU in spanwise direction 
!     and accumulation
      CALL red_z (Wa,red1,1,M1,1,M2,1,M3)  !(U x rotU)_x - Average
      CALL red_z (Wb,red2,1,M1,1,M2,1,M3)  !(U x rotU)_y - Average
      CALL red_z (Wc,red3,1,M1,1,M2,1,M3)  !(U x rotU)_z - Average
      DO i=1,M1
      DO j=1,M2
      cx_acc(i,j)=(cx_acc(i,j)*(No_Sam-1)+red1(i,j))/No_Sam
      cy_acc(i,j)=(cy_acc(i,j)*(No_Sam-1)+red2(i,j))/No_Sam
      cz_acc(i,j)=(cz_acc(i,j)*(No_Sam-1)+red3(i,j))/No_Sam
      ENDDO
      ENDDO

!     dump statistics at last time step (see SUB dump)
!     increase the sample counter
      No_Sam=No_Sam+1

      RETURN
      END   !statistics



!----------------------------------------------------------------------

      SUBROUTINE red_z (W_3d,red_out,imin,imax,jmin,jmax,zlo,zhi)
!     -----------------------------------------------------------
      INCLUDE 'jet3d.inc'

      INTEGER imin,imax,jmin,jmax,zlo,zhi,i,j,k
      REAL W_3d(0:N1,0:N2,0:N3),red_out(imin:imax,jmin:jmax), sum
!     Wd - work array defined in jet3d.inc
!     dz(k) - is defined in jet3d.inc

      DO i=imin,imax
            DO j=jmin,jmax
                  sum=0.0
                  DO k=zlo,zhi
                  sum=sum+W_3d(i,j,k)*dz(k)
                  ENDDO !k
            red_out(i,j)=sum/(z(zhi+1)-z(zlo))
            ENDDO
      ENDDO

      RETURN
      END   !red_z

!----------------------------------------------------------------------

       SUBROUTINE Spectra
!      ------------------
       
    INCLUDE 'jet3d.inc'
    
    INTEGER :: i,j,k
    
!Writes velocity components to files for manipulation
!to produce energy spectra

!Streamwise corellation
!DO I=1,M1
!WRITE(71,*) i, (U1(i,56,32) - u1_acc(i,56)), (U2(i,56,32)-u2_acc(i,56)),&
!            (U3(i,56,32)-u3_acc(i,56)),(T1(i,56,32)-T1_acc(i,56))
!ENDDO

!Cross stream corellation at i=220
!DO J= 10, 90
!WRITE(72,*) j, (U1(220,j,32) - U1_acc(220,j)), (U2(220,j,32)-u2_acc(220,j)), &
!            (U3(220,j,32)-u3_acc(220,j)), (T1(220,j,32)-T1_acc(220,j))
!ENDDO

!Spanwise corellation at i=220, j=56
!DO K=1, M3
!WRITE(73,*) k, (U1(220,56,k)- u1_acc(220,56)), (U2(220,56,k)-u2_acc(220,56)), &
!            (U3(220,56,k)-u2_acc(220,56)), (T1( 220,56,k)-T1_acc(220,56))
!ENDDO

RETURN

END SUBROUTINE Spectra

!----------------------------------------------------------------------


       SUBROUTINE Acoustics
!      --------------------
       INCLUDE 'jet3d.inc'

       INTEGER i,j,k
       REAL vol_int11, vol_int22, vol_int33, vol_int12, vol_int13,     &
            vol_int23, surf_int12, surf_int22, surf_int32 
       REAL u_cc,v_cc,w_cc,du_dy,dv_dy,dw_dy

!      1. integration of volume sources
!      a) initialization
       vol_int11=0.0
       vol_int22=0.0
       vol_int33=0.0
       vol_int12=0.0
       vol_int13=0.0
       vol_int23=0.0
!     
       DO i=1,M1
       DO j=1,M2
       DO k=1,M3
!           b) bring the volicties to the cell center
            u_cc=0.5*(U1(i,j,k)+U1(i+1,j,k))
            v_cc=0.5*(U2(i,j,k)+U2(i,j+1,k))
            w_cc=0.5*(U3(i,j,k)+U3(i,j,k+1))
!           c) integration
            vol_int11=vol_int11+u_cc*u_cc*dx(i)*dy(j)*dz(k)
            vol_int22=vol_int22+v_cc*v_cc*dx(i)*dy(j)*dz(k)
            vol_int33=vol_int33+w_cc*w_cc*dx(i)*dy(j)*dz(k)
            vol_int12=vol_int12+u_cc*v_cc*dx(i)*dy(j)*dz(k) 
            vol_int13=vol_int13+u_cc*w_cc*dx(i)*dy(j)*dz(k)
            vol_int23=vol_int23+v_cc*w_cc*dx(i)*dy(j)*dz(k) 
       ENDDO
       ENDDO
       ENDDO
!
!      2. integration of surface sources
!      a) initialization
       surf_int12=0.0
       surf_int22=0.0
       surf_int32=0.0
!
       j=1  !surface
       DO i=1,M1 
       DO k=1,M3
!            b) compute the wall normal velocity gradients at wall,
!               gradients defined at the center of the cell bottom (i,k,j=1)
             du_dy=0.5*(U1(i+1,j,k)-U1(i+1,j-1,k)+U1(i,j,k)-           &
                   U1(i,j-1,k))*dy1(j)
             dw_dy=0.5*(U3(i,j,k+1)-U3(i,j-1,k+1)+U3(i,j,k)-           &
                   U3(i,j-1,k))*dy1(j)
             dv_dy=0.5*(U2(i,j+1,k)-U2(i,j,k))*dy0(j) !NOTE grad(j=0,1)=0.0
!                                                      & dy0(0)=dy0(1)
!            c) integration
             surf_int12=surf_int12+du_dy*dx(i)*dz(k)
             surf_int22=surf_int22+dv_dy*dx(i)*dz(k)
             surf_int32=surf_int32+dw_dy*dx(i)*dz(k)
       ENDDO
       ENDDO
!
!      3. write results to the file
       WRITE (12,100) Time, vol_int11, vol_int22, vol_int33, vol_int12,&
              vol_int13, vol_int23, surf_int12, surf_int22, surf_int32

100    FORMAT (e14.6,e18.10,e18.10,e18.10,e18.10,e18.10,e18.10,e18.10, &
       e18.10,e18.10)
       RETURN
       END  !Acoustics
!----------------------------------------------------------------------

!       SUBROUTINE Trace
!      ----------------
!       INCLUDE 'jet3d.inc'
!       INTEGER i,lognum
!       REAL u_cc,v_cc,w_cc

!       DO i=1,Max_Traces
!       lognum=30+i
!      1. bring the velocities to the cell center
!       u_cc=0.5*(U1(Xtr(i),Ytr(i),Ztr(i))+U1(Xtr(i)+1,Ytr(i),Ztr(i)))
!       v_cc=0.5*(U2(Xtr(i),Ytr(i),Ztr(i))+U2(Xtr(i),Ytr(i)+1,Ztr(i)))
!       w_cc=0.5*(U3(Xtr(i),Ytr(i),Ztr(i))+U3(Xtr(i),Ytr(i),Ztr(i)+1))
!      2. write them to file
!       WRITE (lognum,105) Time,u_cc,v_cc,w_cc,P1(Xtr(i),Ytr(i),Ztr(i))
!       ENDDO
!
!105    FORMAT (e14.6,e18.10,e18.10,e18.10,e18.10)
!       RETURN
!       END  !Trace
!----------------------------------------------------------------------
 
       SUBROUTINE Finish
!      -----------------
         INCLUDE 'jet3d.inc'

         PRINT *, Abort_Message, ' at step number ', Step_No
!        CALL Dump
         CALL Shut_Files

         RETURN
       END

!----------------------------------------------------------------------

      SUBROUTINE Assign_Files
!     -----------------------
      INCLUDE 'jet3d.inc'
      INTEGER i,lognum
      CHARACTER*1 ch1
      CHARACTER*12 tracefile

!     open all files for monitoring and accumulation
      OPEN (UNIT=9,FILE='monitor.3d')
      REWIND (9)     
      WRITE (9,*) '------------------ M O N I T O R -------------------'
      WRITE (9,*) 'step  time  CFLmax Cox Coy Coz i j k INTdiv(U)dV  CPU- total/pressure'
!     *********************** end monitor
      OPEN (UNIT=12,FILE='sound.3d')
      REWIND (12)     
      WRITE (12,*) '---------------- A C O U S T I C S ----------------'
      WRITE (12,*) 'time, v_int. (uu,vv,ww,uv,uw,vw), s_int. (uy,vy,wy)' 
!     *********************** end acoustics
!      DO i=1,Max_Traces
!      ch1=CHAR(i+48)      !ch1=1...9 only
!      lognum=30+i
!      tracefile='trace'//ch1//'.3d'
!      OPEN (UNIT=lognum,FILE=tracefile)
!      REWIND (lognum)
!      WRITE (lognum,*) '------------ T R A C E P O I N T --------------'
!      WRITE (lognum,*) 'time,       U1          U2          U3       P1'
!      ENDDO
!     *********************** end trace points

      RETURN
      END      !Assign_Files

!----------------------------------------------------------------------

         SUBROUTINE Shut_Files
!        ---------------------
         INCLUDE 'jet3d.inc'
         INTEGER i,lognum

!        close all files for monitoring and accumulation  
         WRITE (9,101) 'END'
         CLOSE (UNIT=9, STATUS='KEEP')    !Monitor

         WRITE (12,101) 'END'
         CLOSE (UNIT=12, STATUS='KEEP')    !Acoustics
         
!         DO i=1,Max_Traces
!         lognum=30+i
!         WRITE (lognum,101) 'END'
!         CLOSE (UNIT=lognum, STATUS='KEEP')     !Trace Points
!         ENDDO

101      FORMAT (a3)   
         RETURN
         END      !Shut_files

!----------------------------------------------------------------------

       SUBROUTINE X_Table ( Channel, delt )
!      ------------------------------------
         INCLUDE 'jet3d.inc'
         INTEGER Channel, i , idum
         REAL delt, xtemp(0:Q1+2)

!         Channel = Channel/1 ! prevents a warning
         IF( Channel .GT. 1 ) GOTO 11
         Delta_x = delt
         Size_x = delt * M1
         Rdx = 1.0 / delt
         DO i = 0, N1
           dx(i) = delt    ! used by pressure solver
           dx1(i) = 1.0/dx(i)
         ENDDO
!PVM
         Offset_x = delt * Q1 / L1 * O1
         x(0) = -delt + Offset_x
         DO 10, i = 0+1, N1+1
           x(i) = x(i-1) + delt
 10      CONTINUE
!         PRINT *, 'X_Mesh is uniform'

 11      CONTINUE
         IF( Channel .NE. 0 ) THEN
           OPEN( Channel,FILE='xcoord.3d')               
                 
!PVM
           DO i = 1, M1+1
           READ( Channel, * ) idum, xtemp(i)
           IF( i .NE. idum )THEN
           WRITE(6,*)' Incompatible x-mesh! '
           ENDIF
           ENDDO
!
           CLOSE (Channel)
!
           xtemp(0)=2.0*xtemp(1)-xtemp(2)
           xtemp(Q1+2)=2.0*xtemp(Q1+1)-xtemp(Q1)
           DO i = 0, N1+1
             x(i) = xtemp(i+O1)
           ENDDO
           DO i = 0, N1
             dx(i) = x(i+1) - x(i)
           ENDDO
           DO i = 1, N1
             dx0(i) = 1.0/dx(i)
             dx1(i) = 2.0/(dx(i) + dx(i-1))
           ENDDO
           dx0(0) = 1.0/dx(0)
           dx1(0) = dx1(1)

          DO I = 0, M1
             xv(i) = 0.5 * (x(i) + x(i+1))
			 xv(N1) = xv(M1) + dx(M1)
          ENDDO

		  DO I=0, M1
			 dxv(i) = (xv(i+1) - xv(i))
		  ENDDO
		  dxv(N1) = dxv(M1)			  
         ENDIF


         RETURN
 
! 30      CONTINUE
         PRINT *, 'End of file on x mesh'
         Abort = .TRUE.

         RETURN
       END
 
       SUBROUTINE Y_Table ( Channel , delt )
!      -------------------------------------
         INCLUDE 'jet3d.inc'
         INTEGER Channel, j , jdum
         REAL delt, ytemp(0:Q2+2)

!         delt = delt / 1 ! prevents a warning
         OPEN( Channel, FILE='ycoord.3d' )
!
!PVM
         DO j=1,m2+1
         READ ( Channel, * ) jdum, ytemp(j)
         IF( j .NE. jdum )THEN
         WRITE(6,*) ' Incompatible y-mesh !'
         STOP
         ENDIF
         ENDDO
!
         CLOSE (Channel)
         ytemp(0)=2*ytemp(1)-ytemp(2)
         ytemp(Q2+2)=2*ytemp(Q2+1)-ytemp(Q2)
         DO j = 0, N2+1
           y(j) = ytemp(j+O2)
         ENDDO
         DO j = 0, N2
           dy(j) = y(j+1) - y(j)
         ENDDO
         DO j = 1, N2
           cy0(j) = dy(j-1) / ( dy(j-1) + dy(j) )
           cy1(j) = dy(j) / ( dy(j-1) + dy(j) )
           dy0(j) = 1.0 / dy(j)
           dy1(j) = 2.0 / ( dy(j-1) + dy(j) )
         ENDDO
         cy0(0) = cy0(1)
         cy1(0) = cy1(1)
         dy0(0) = 1.0 / dy(0)
         dy1(0) = dy1(1)

         DO J=0, N2
            yu(j) = 0.5 * (y(j) + y(j+1))
            yu(N2) = yu(M2) + dy(M2)
         ENDDO

		DO J=0, M2
		     dyu(j) = (yu(j+1) - yu(j))
		ENDDO
		dyu(N2) = dyu(M2)
         
         RETURN
 
! 20     CONTINUE
         Abort_Message = 'End of file on y mesh.'
         PRINT *, '*** ', Abort_Message
         Abort = .TRUE.

         RETURN
       END

       SUBROUTINE Z_Table ( Channel, delt )
!      ------------------------------------
         INCLUDE 'jet3d.inc'
         INTEGER Channel, k
         REAL delt

!         Channel = Channel/1 ! prevents a warning
         Delta_z = delt
         Size_z = delt * M3
         Rdz = 1.0 / delt

         DO k = 0, N3
           dz(k) = delt
         ENDDO
!PVM
         Offset_z = delt * Q3 / L3 * O3
         z(0) = -delt + Offset_z
         DO k = 1, N3+1
           z(k) = z(k-1) + delt
         ENDDO
         

         RETURN

       END

!----------------------------------------------------------------------
!----------------------------------------------------------------------
 
      SUBROUTINE fft_Forward(A)
!     -------------------------
!     Bit-scrambled radix-2 fft of 3D real array A, in z dimension.
!     Uses odd j as real part, even j as imaginary part.
 
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'
      REAL A(0:N1,0:N2,0:N3), B(1:M1,1:M2)
      REAL Rs,Is,real,imag
      REAL real1,imag1,real2,imag2 !: used for Herm wrap only
      INTEGER k1, k2 !: used for Herm wrap only   
      INTEGER max, s, pp, c, m, i, j, k, kk, l, jump
      LOGICAL upper

      IF (Q3 .LE. 1) RETURN     

!PVM  Layers down to base-8, extensive mods for PVM
      IF (Q3 .GT. 4) THEN
      max = Q3/2
 10   CONTINUE ! max loop
        s = 2*max
        pp = Q3/s
      IF ( max .LT. M3 ) THEN  
!       Standard fft algorithm
        c = 0
        DO m = 1, max
          real = Re(c)
          imag = Im(c)
          c = c+pp
          DO k = m, M3, s
            kk = k+max
          DO i = 1, M1
          DO j = 1, M2, 2
            l = j+1
            Rs = (A(i,j,k)-A(i,j,kk))/2
            Is = (A(i,l,k)-A(i,l,kk))/2
            A(i,j,k)  = A(i,j,k)-Rs
            A(i,l,k)  = A(i,l,k)-Is
            A(i,j,kk) = real*Rs+imag*Is
            A(i,l,kk) = real*Is-imag*Rs
          ENDDO ! j
          ENDDO ! i
          ENDDO ! k
        ENDDO ! m
      ELSE      ! max .GE. M3
!PVM    Data exchange with other processors
        jump = max/M3
        upper = (MOD(me(3)/jump,2) .NE. 0) 
        c = pp * M3 * MOD(me(3),jump)
        DO k = 1, M3
          real = Re(c)
          imag = Im(c)
          c = c+pp
!          IF (upper) THEN
!           talk with lower processor
!            msgtype = me(1)*1000+me(2)*100-jump
!            CALL pvmfinitsend(PVMDEFAULT, info)
!            DO j = 1, M2
!              CALL pvmfpack( REAL8, A(1,j,k), M1, 1, info )
!            ENDDO
!            CALL pvmfsend( Ztid(-jump), msgtype, info ) 
!            msgtype = me(1)*1000+me(2)*100+jump
!            CALL pvmfrecv( Ztid(-jump), msgtype, info ) 
!            DO j = 1, M2
!              CALL pvmfunpack( REAL8, B(1,j), M1, 1, info )
!            ENDDO
            DO i = 1, M1
            DO j = 1, M2, 2
              l = j+1
              Rs = (B(i,j)-A(i,j,k))/2
              Is = (B(i,l)-A(i,l,k))/2
              A(i,j,k) = real*Rs+imag*Is
              A(i,l,k) = real*Is-imag*Rs
            ENDDO 
            ENDDO 
!          ELSE ! lower
!           talk with upper processor
!            msgtype = me(1)*1000+me(2)*100+jump
!            CALL pvmfinitsend(PVMDEFAULT, info)
!            DO j = 1, M2
!              CALL pvmfpack( REAL8, A(1,j,k), M1, 1, info )
!            ENDDO
!            CALL pvmfsend( Ztid(jump), msgtype, info )
!            msgtype = me(1)*1000+me(2)*100-jump
!            CALL pvmfrecv( Ztid(jump), msgtype, info )
!            DO j = 1, M2
!              CALL pvmfunpack( REAL8, B(1,j), M1, 1, info )
!            ENDDO
            DO i = 1, M1
            DO j = 1, M2 !, 2
!              l = j+1
              A(i,j,k) = (A(i,j,k)+B(i,j))/2
!              A(i,l,k) = (A(i,l,k)+B(i,l))/2  presumably redundant
            ENDDO 
            ENDDO 
!          ENDIF ! upper else lower
        ENDDO ! k (= m)
      
      ENDIF ! max .LT. M3

        max = max/2
      IF ( max .GT. 2) GOTO 10  ! end of max loop
      ENDIF
 
!     base-4 layer
      IF (Q3 .GT. 2) THEN
      DO k = 1, M3, 4
        kk = k+2
      DO i = 1, M1
      DO j = 1, M2, 2
        l = j+1
        A(i,j,kk) = (A(i,j,k)-A(i,j,kk))/2
        A(i,l,kk) = (A(i,l,k)-A(i,l,kk))/2
        A(i,j,k)  =  A(i,j,k)-A(i,j,kk)
        A(i,l,k)  =  A(i,l,k)-A(i,l,kk)
      ENDDO
      ENDDO
      ENDDO
      DO k = 2, M3, 4
        kk = k+2
      DO i = 1, M1
      DO j = 1, M2, 2
        l = j+1
        Rs  = (A(i,j,k)-A(i,j,kk))/2
        A(i,j,kk) = (A(i,l,kk)-A(i,l,k))/2
        A(i,j,k)  = A(i,j,k)-Rs
        A(i,l,k)  = A(i,l,k)+A(i,j,kk)
        A(i,l,kk) = Rs
      ENDDO
      ENDDO
      ENDDO
      ENDIF
 
!     base-2 layer
      DO k = 1, M3, 2
        kk = k+1
      DO i = 1, M1
      DO j = 1, M2, 2
        l = j+1
        A(i,j,kk) = (A(i,j,k)-A(i,j,kk))/2
        A(i,l,kk) = (A(i,l,k)-A(i,l,kk))/2
        A(i,j,k)  =  A(i,j,k)-A(i,j,kk)
        A(i,l,k)  =  A(i,l,k)-A(i,l,kk)
      ENDDO
      ENDDO
      ENDDO
 
!     Hermitian/skew-Hermitian unwrap
      DO k = 2, M3/2
!       look up current (bit-reversed) positions of +/- k frequency
        k1 = posf(k)
        k2 = negf(k)
!       put Hermitian transform in A(j) and 
!       skew-Hermitian transform in A(l).
        DO i = 1, M1
        DO j = 1, M2, 2
          l = j+1
          real1 = A(i,j,k1)
          imag1 = A(i,l,k1)
          real2 = A(i,j,k2)
          imag2 = A(i,l,k2)
          A(i,j,k1) = ( real1+real2 )/2  ! real part of A(j)
          A(i,j,k2) = ( imag1-imag2 )/2  ! imag part of A(j)
          A(i,l,k1) = ( imag1+imag2 )/2  ! real part of A(l)
          A(i,l,k2) = ( real2-real1 )/2  ! imag part of A(l)
        ENDDO
        ENDDO
      ENDDO
 
      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE fft_Inverse(A)
!     ------------------------
!     Bit-scrambled radix-2 inverse z-fft of 3D Hermitian sequences A.
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'
      REAL A(0:N1,0:N2,0:N3), B(0:N1,0:N2)
      REAL Rt,It,real,imag
      REAL real1,imag1,real2,imag2 !: used for Herm wrap only   
      INTEGER k1, k2 !: used for Herm wrap only   
      INTEGER max, s, pp, c, m, i, j, k, kk, l, jump
      LOGICAL upper

      IF (Q3 .LE. 1) RETURN  

!     Hermitian/skew-Hermitian wrap
      DO k = 2, M3/2
        k1 = posf(k)
        k2 = negf(k)
        DO i = 1, M1
        DO j = 1, M2, 2
          l = j+1
          real1 = A(i,j,k1)
          imag1 = A(i,l,k1)
          real2 = A(i,j,k2)
          imag2 = A(i,l,k2)
          A(i,j,k1) = real1-imag2
          A(i,j,k2) = real1+imag2
          A(i,l,k1) = imag1+real2
          A(i,l,k2) = imag1-real2
        ENDDO
        ENDDO
      ENDDO
 
!     base-2 layer fft: should be local
      DO k = 1, M3, 2
        kk = k+1
      DO i = 1, M1
      DO j = 1, M2, 2
        l = j+1
        Rt = A(i,j,kk)
        It = A(i,l,kk)
        A(i,j,kk) = A(i,j,k)-Rt
        A(i,l,kk) = A(i,l,k)-It
        A(i,j,k)  = A(i,j,k)+Rt
        A(i,l,k)  = A(i,l,k)+It
      ENDDO
      ENDDO
      ENDDO
 
!     base-4 layer fft: should be local
      IF (Q3 .GT. 2) THEN
      DO k = 1, M3, 4
        kk = k+2
      DO i = 1, M1
      DO j = 1, M2, 2
        l = j+1
        Rt = A(i,j,kk)
        It = A(i,l,kk)
        A(i,j,kk) = A(i,j,k)-Rt
        A(i,l,kk) = A(i,l,k)-It
        A(i,j,k)  = A(i,j,k)+Rt
        A(i,l,k)  = A(i,l,k)+It
      ENDDO
      ENDDO
      ENDDO
      DO k = 2, M3, 4
      kk = k+2
      DO i = 1, M1
      DO j = 1, M2, 2
        l = j+1
        Rt = A(i,l,kk)
        It = A(i,j,kk)
        A(i,j,kk) = A(i,j,k)-Rt
        A(i,l,kk) = A(i,l,k)+It
        A(i,j,k)  = A(i,j,k)+Rt
        A(i,l,k)  = A(i,l,k)-It
      ENDDO
      ENDDO
      ENDDO
      ENDIF
 
!PVM  base-8 and upwards fft, extensive mods for PVM
      IF (Q3 .GT. 4) THEN
      max=4
 10   CONTINUE ! max loop
        s = 2*max
        pp = Q3/s
      IF ( max .LT. M3 ) THEN  
!       Standard fft algorithm
        c = 0
        DO m = 1, max
          real = Re(c)
          imag = Im(c)
          c = c+pp
          DO k = m, M3, s
          kk = k+max
          DO i = 1, M1
          DO j = 1, M2, 2
            l = j+1
            Rt = real*A(i,j,kk)-imag*A(i,l,kk)
            It = real*A(i,l,kk)+imag*A(i,j,kk)
            A(i,j,kk) = A(i,j,k)-Rt
            A(i,l,kk) = A(i,l,k)-It
            A(i,j,k)  = A(i,j,k)+Rt
            A(i,l,k)  = A(i,l,k)+It
          ENDDO ! j
          ENDDO ! i
          ENDDO ! k
        ENDDO   ! m
      ELSE      ! max .GE. M3
!PVM    Data exchange with other processors
        jump = max/M3
        upper = (MOD(me(3)/jump,2) .NE. 0) 
        c = pp * M3 * MOD(me(3),jump)
        DO k = 1, M3
          real = Re(c)
          imag = Im(c)
          c = c+pp
!          IF (upper) THEN
!           talk with lower processor
!            msgtype = me(1)*1000+me(2)*100-jump
!            CALL pvmfinitsend(PVMDEFAULT, info)
!            DO j = 1, M2
!              CALL pvmfpack( REAL8, A(1,j,k), M1, 1, info )
!            ENDDO
!            CALL pvmfsend( Ztid(-jump), msgtype, info ) 
!            msgtype = me(1)*1000+me(2)*100+jump
!            CALL pvmfrecv( Ztid(-jump), msgtype, info ) 
!            DO j = 1, M2
!              CALL pvmfunpack( REAL8, B(1,j), M1, 1, info )
!            ENDDO
            DO i = 1, M1
            DO j = 1, M2, 2
              l = j+1
              A(i,j,k) = B(i,j)-(real*A(i,j,k)-imag*A(i,l,k)) 
              A(i,l,k) = B(i,l)-(real*A(i,l,k)+imag*A(i,j,k))
            ENDDO ! j
            ENDDO ! i
!          ELSE ! lower
!           talk with upper processor
!            msgtype = me(1)*1000+me(2)*100+jump
!            CALL pvmfinitsend(PVMDEFAULT, info)
!            DO j = 1, M2
!              CALL pvmfpack( REAL8, A(1,j,k), M1, 1, info )
!            ENDDO
!            CALL pvmfsend( Ztid(+jump), msgtype, info ) 
!            msgtype = me(1)*1000+me(2)*100-jump
!            CALL pvmfrecv( Ztid(+jump), msgtype, info ) 
!            DO j = 1, M2
!              CALL pvmfunpack( REAL8, B(1,j), M1, 1, info )
!            ENDDO
            DO i = 1, M1
            DO j = 1, M2, 2
              l = j+1
              A(i,j,k) = A(i,j,k)+real*B(i,j)-imag*B(i,l)
              A(i,l,k) = A(i,l,k)+real*B(i,l)+imag*B(i,j)
            ENDDO ! j
            ENDDO ! i
!          ENDIF ! upper / lower
        ENDDO ! k
      ENDIF ! max .LT. M3
        max = s
      IF (Q3 .GT. max) GOTO 10
      ENDIF
 
      RETURN
      END

!----------------------------------------------------------------------
 
      SUBROUTINE fft_Table
!     --------------------
!     (pseudo-wave-numbers)**2 = ksqu, and frequency lookup tables

      INCLUDE 'jet3d.inc'
      REAL Pi 
      PARAMETER (Pi = 3.14159265308979)
      INTEGER j, k, m, q, wave
      REAL phi, temp

!PVM  Q3
      DO k = 1, Q3
        wave = k-1
        IF (k .GT. M3/2+1) wave = M3-wave
        posf(k) = k
        negf(k) = M3-k+2
        IF (k .EQ. 1) negf(k)=1
        ksqu(k) = (2.0*SIN(Pi*wave/M3))**2/(Delta_z*Delta_z)
      ENDDO
!     ...in bit-reversed order.
      j = 1
      DO k = 1, Q3
        IF (j .GT. k) THEN
          temp = ksqu(k)
          ksqu(k) = ksqu(j)
          ksqu(j) = temp
          DO q = 1, M3
           IF (posf(q) .EQ. k ) THEN
             posf(q)=j
           ELSE
             IF (posf(q) .EQ. j) posf(q)=k
           ENDIF
           IF (negf(q) .EQ. k) THEN
             negf(q)=j
           ELSE
             IF (negf(q) .EQ. j) negf(q)=k
           ENDIF
          ENDDO
        ENDIF
        m = Q3/2
 10     CONTINUE
          IF (m .LT. 1 .OR. j .LE. m) GOTO 20
          j = j-m
          m = m/2
        GOTO 10
 20     CONTINUE
        j = j+m
      ENDDO ! k
 
!     trig factors for fft.
      DO k = 1, Q3/2
        phi = (1-k)*2.0*Pi/M3
        Re(k-1) = COS(phi)
        Im(k-1) = SIN(phi)
      ENDDO
 
      RETURN
      END
!-------------------------------------------------------------

      SUBROUTINE Ini_U_P
!     ------------------
      INCLUDE 'jet3d.inc'
      INTEGER i,j,k

        DO k = 0, N3
        DO j = 0, N2
        DO i = 0, N1
          U1(i,j,k) = 0.0
          U2(i,j,k) = 0.0
          U3(i,j,k) = 0.0
          P1(i,j,k) = 0.0
          T1(i,j,k) = 0.0
          pass1(i,j,k) = 0.0
          pasf1(i,j,k) = 0.0
        ENDDO
        ENDDO
        ENDDO

      RETURN
      END

      SUBROUTINE Ini_W
!     ----------------
      INCLUDE 'jet3d.inc'
      INTEGER i,j,k

        DO k = 0, N3
        DO j = 0, N2
        DO i = 0, N1
          Wa(i,j,k) = 0.0
          Wb(i,j,k) = 0.0
          Wc(i,j,k) = 0.0
          Wd(i,j,k) = 0.0
          We(i,j,k) = 0.0
          Wg(i,j,k) = 0.0
          Wh(i,j,k) = 0.0
        ENDDO
        ENDDO
        ENDDO

      RETURN
      END

      SUBROUTINE Ini_H
!     ----------------
      INCLUDE 'jet3d.inc'
      INTEGER i,j,k

        DO k = 0, N3
        DO j = 0, N2
        DO i = 0, N1
          H1(i,j,k) = 0.0
          H2(i,j,k) = 0.0
          H3(i,j,k) = 0.0
          R1(i,j,k) = 0.0
        ENDDO
        ENDDO
        ENDDO

      RETURN
      END

      SUBROUTINE Ini_S
!     ----------------
      INCLUDE 'jet3d.inc'
      INTEGER i,j

        DO j = 1, M2
        DO i = 1, M1
          u1_acc(i,j) = 0.0
          u2_acc(i,j) = 0.0
          u3_acc(i,j) = 0.0
          p1_acc(i,j) = 0.0
          t1_acc(i,j) = 0.0
          uu_acc(i,j) = 0.0
          vv_acc(i,j) = 0.0
          ww_acc(i,j) = 0.0
          uv_acc(i,j) = 0.0
          uw_acc(i,j) = 0.0
          vw_acc(i,j) = 0.0
          up_acc(i,j) = 0.0
          vp_acc(i,j) = 0.0
          wp_acc(i,j) = 0.0
          cx_acc(i,j) = 0.0
          cy_acc(i,j) = 0.0
          cz_acc(i,j) = 0.0
        ENDDO
        ENDDO

      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE int2char (inum,charout)
!     ----------------------------------
!     converts the INTEGER number inum of ndig digits into a character
!     string, max. length of the character*1*(lmax)
!     uses power_10, an integer-function
!
      PARAMETER (lmax=6)
      LOGICAL chflag
      CHARACTER*1 charout(lmax)
      INTEGER basis(lmax+1),qi(lmax),ndig,inum,ii,ll
      INTEGER prodi,asci,power_10
      REAL qr
!
!     calculation of the length of the integer (number of digits)
      ll=0
      qr=inum/power_10(ll)
      IF (qr .LE. 1.0) THEN 
            ndig=1
      ELSE
            ll=1
55          qr=inum/power_10(ll)
            IF (qr .LT. 1.0) THEN 
                  ndig=ll
                  chflag=.FALSE.
            ELSE
                  ll=ll+1
                  chflag=.TRUE.
            ENDIF
            IF (chflag) GOTO 55
      ENDIF
!
!     separation of digits
      basis(1)=inum
      DO 111 ii=1,ndig
            ll=ndig-ii
            qr=basis(ii)/power_10(ll)
            qi(ii)=INT(qr)
            prodi=qi(ii)*power_10(ll)
            basis(ii+1)=basis(ii)-prodi
111   CONTINUE
!
!     initialization of string
      DO 33 ii=1,lmax
            charout(ii)='_'
33    CONTINUE
!
!     convert digits into chararcter codes
      DO 222 ii=1,ndig
            asci=qi(ii)+48
            charout(ii)=char(asci)
222   CONTINUE
!
      RETURN
      END
!     ################### end of int2char #################
!
      INTEGER FUNCTION power_10 (jj)
!     ------------------------------
!     calculates 10 to the jj power
!
      INTEGER jj,kk
!
      power_10=1
      DO 333 kk=1,jj
            power_10=power_10*10
333   CONTINUE
!
      RETURN
      END
!     ################### end of power_10 #################
!
!     PVM - Routines
!     --------------
!
      SUBROUTINE pvmfpack (int0,real1,int1,int2,int3)
!     -----------------------------------------------
!
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'
!
      REAL real1
      INTEGER int0,int1,int2,int3
!
      RETURN
      END
!     ################### end of pvmfpack #################
!
      SUBROUTINE pvmfsend (int0,int1,int2)
!     ------------------------------------
!
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'    
!
      INTEGER  int0,int1,int2
!
      RETURN
      END
!     ################### end of pvmfsend #################
!     
      SUBROUTINE pvmfunpack (int0,real1,int1,int2,int3)
!     -------------------------------------------------
!
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'    
!
      INTEGER int0,int1,int2,int3
      REAL real1
!
      RETURN
      END
!     ################## end of pvmfunpack ################
!     
      SUBROUTINE pvmfinitsend (int0,int1)
!     -----------------------------------
!
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'    
!
      INTEGER int0,int1
!
      RETURN
      END
!     ################ end of pvmfinitsend ################
!     
      SUBROUTINE pvmfrecv (int0,int1,int3)
!     ------------------------------------
!
      INCLUDE 'jet3d.inc'
      INCLUDE 'pvm.inc'    
!
      INTEGER int0,int1,int3
!
      RETURN
      END
!     ################## end of pvmfrecv ##################
