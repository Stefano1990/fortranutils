       SUBROUTINE dataread 
!      -------------------- 
                   
         implicit none
            
!include files here
#include "dim.inc"
#include "cb.control.inc"
#include "cb.comm.inc"

         character*30 string30, temp30
         character* 6 string, temp
         real value
         integer round, i, j, k

         if(master) then 
         write (*,*) 'DATA_READ subroutine...' 
         endif
         OPEN ( 10, FILE='input.dat', STATUS='OLD' )

 10      continue
           read (10, '(F20.10, A)', end=20 ) value, string30
           if ( ABS(Value) .gt. 0.0 ) then
           if(master) then
             print *, string30, '  =  ', Value
!           else
!             print *, string30
           endif
           endif
           Round = INT( Value+0.000001 )
           
           if (string30 .eq. ' ') string30 = '*'
           do while ( string30(1:1) .eq. ' ' )
             Temp30 = string30(2:30) // ' '
             string30 = Temp30
           enddo
           i = INDEX( string30,' ' )
           if (i .gt. 0) then 
             Temp30 = string30(1:i-1)
             string30 = Temp30
           endif
           string = string30(1:6)
           do i = 1, 6
             k = ICHAR(string(i:i))
             if (k .gt. 95) then
               k = k-32
               Temp = string(1:i-1) // CHAR(k) // string(i+1:6)
               string = Temp
             endif
           enddo

             if ( string .eq. 'VISCOS') then
               visc = value 
             ENDIF
             if ( string .eq. 'LAST' )  last = Round
             if ( string .eq. 'TIME' )  dtime = Value
             if ( string .eq. 'NEW' )   New_Run = .TRUE.
             if ( string .eq. 'OLD' )   New_Run = .FALSE.
             if ( string .eq. 'DEBUG' ) debugging = .true.
             if ( string .eq. 'DUMP' )  Dump_Freq = Round
             if ( string .eq. 'SPLITT' )  jsplit = Round
             if ( string .eq. 'DIRECT') then
                LES = .false.
                Model = None
             endif
             if ( string .eq. 'SMAGOR') then
               LES = .true.
               Model = Smagorinsky
               csles = Value     
             endif
              if ( string .eq. 'WALE') then
               LES = .true.
               Model = Wale
               csles = Value     
             endif            
             if ( string .eq. 'RHOHI' ) rhohi = Value
             if ( string .eq. 'RHOLOW') rholo = Value
             if ( string .eq. 'NUMVAR') numvar = Round
             if ( string .eq. 'TRELAX')  t_relax = Value
             if ( string .eq. 'ROTATE') rotate = Round
             if ( string .eq. 'U-CYC' ) ncycle(1) = Round
             if ( string .eq. 'V-CYC' ) ncycle(2) = Round
             if ( string .eq. 'W-CYC' ) ncycle(3) = Round
             if ( string .eq. 'P-CYC' ) ncycle(4) = Round
             if ( string .eq. 'SC-CYC') ncycle(5) = Round 
             if ( string .eq. 'NUMSOL') numsolve = Round 
             if ( string .eq. 'RELU'  ) relax(1) = Round 
             if ( string .eq. 'RELV'  ) relax(2) = Round
             if ( string .eq. 'RELW'  ) relax(3) = Round
             if ( string .eq. 'RELP'  ) relax(4) = Round
             if ( string .eq. 'RELSC'  ) relax(5) = Round
             if ( string .eq. 'PRESBL' ) refpresblock = Round
             if ( string .eq. 'INLET1' ) inlet1 = Round
             if ( string .eq. 'INLET2' ) inlet2 = Round
             if ( string .eq. 'INLET3' ) inlet3 = Round
             if ( string .eq. 'INLET4' ) inlet4 = Round
             if ( string .eq. 'RHOINI' ) rhoini = Value
             if ( string .eq. 'WRITEW' ) writewall = .true.
                           
             
! Statistics and flow viz outputs
             if ( String .EQ. 'STATIS') Stat_Freq = Round
             if ( String .EQ. 'FIELDS') Fiel_Freq = Round
             if ( String .EQ. 'SPECTR') Spec_Freq = Round
             if ( String .EQ. 'STATST') statstart = Round
             if ( String .EQ. 'ADDBAC') addback = Round
!output in original ifile format - temporary
!as will change to single file in the future			 
             if ( String .EQ. 'U' )      ifile(1) = .TRUE.
             if ( String .EQ. 'V' )      ifile(2) = .TRUE.
             if ( String .EQ. 'W' )      ifile(3) = .TRUE.
             if ( String .EQ. 'P' )      ifile(4) = .TRUE.
             if ( String .EQ. 'SC1' )    ifile(5) = .TRUE.
             if ( String .EQ. 'HU' )     ifile(11) = .TRUE.
             if ( String .EQ. 'HV' )     ifile(12) = .TRUE.
             if ( String .EQ. 'HW' )     ifile(13) = .TRUE.
             if ( String .EQ. 'HSC1' )   ifile(14) = .TRUE.
             if ( String .EQ. 'RHO' )    ifile(37) = .TRUE.
             if ( String .EQ. 'GAM' )    ifile(42) = .TRUE. 

             if ( string .eq. 'YZ1') yz(1) = Value
             if ( string .eq. 'YZ2') yz(2) = Value
             if ( string .eq. 'YZ3') yz(3) = Value
             if ( string .eq. 'YZ4') yz(4) = Value
             if ( string .eq. 'YZ5') yz(5) = Value
             if ( string .eq. 'YZ6') yz(6) = Value
             if ( string .eq. 'YZ7') yz(7) = Value
             if ( string .eq. 'YZ8') yz(8) = Value
             if ( string .eq. 'YZ9') yz(9) = Value


         if ( string .ne. 'END' ) GOTO 10

         if(master) then 
         WRITE (*,*) '... DATA_READ completed!'
         endif
	 RETURN

 20      PRINT *, '** Warning: abnormal termination of input file.'
         PRINT *, '** END missing from input file.'

         CLOSE (10)



         RETURN
 
       END


