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
             if ( string .eq. 'SIGMA' ) then
              LES = .true.
              Model = Sigma
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
             if ( String .EQ. 'FLDSTR') fieldstart = Round
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
             if ( string .eq. 'YZ10') yz(10) = Value
             if ( string .eq. 'YZ11') yz(11) = Value
             if ( string .eq. 'YZ12') yz(12) = Value
             if ( string .eq. 'YZ13') yz(13) = Value
             if ( string .eq. 'YZ14') yz(14) = Value
             if ( string .eq. 'YZ15') yz(15) = Value
             if ( string .eq. 'YZ16') yz(16) = Value
             if ( string .eq. 'YZ17') yz(17) = Value
             if ( string .eq. 'YZ18') yz(18) = Value
             if ( string .eq. 'YZ19') yz(19) = Value
             if ( string .eq. 'YZ20') yz(20) = Value

             if ( string .eq. 'YPLANE' ) yplane = Value
             if ( string .eq. 'X01') xo(1) = Value
             if ( string .eq. 'X02') xo(2) = Value
             if ( string .eq. 'X03') xo(3) = Value
             if ( string .eq. 'X04') xo(4) = Value
             if ( string .eq. 'X05') xo(5) = Value
             if ( string .eq. 'X06') xo(6) = Value
             if ( string .eq. 'X07') xo(7) = Value
             if ( string .eq. 'X08') xo(8) = Value
             if ( string .eq. 'X09') xo(9) = Value 
             if ( string .eq. 'X10') xo(10) = Value
             if ( string .eq. 'X11') xo(11) = Value
             if ( string .eq. 'X12') xo(12) = Value
             if ( string .eq. 'X13') xo(13) = Value
             if ( string .eq. 'X14') xo(14) = Value
             if ( string .eq. 'X15') xo(15) = Value
             if ( string .eq. 'X16') xo(16) = Value
             if ( string .eq. 'X17') xo(17) = Value
             if ( string .eq. 'X18') xo(18) = Value
             if ( string .eq. 'X19') xo(19) = Value
             if ( string .eq. 'X20') xo(20) = Value
             
             if ( string .eq. 'Y01') yo(1) = Value
             if ( string .eq. 'Y02') yo(2) = Value
             if ( string .eq. 'Y03') yo(3) = Value
             if ( string .eq. 'Y04') yo(4) = Value
             if ( string .eq. 'Y05') yo(5) = Value
             if ( string .eq. 'Y06') yo(6) = Value
             if ( string .eq. 'Y07') yo(7) = Value
             if ( string .eq. 'Y08') yo(8) = Value
             if ( string .eq. 'Y09') yo(9) = Value
             if ( string .eq. 'Y10') yo(10) = Value
             if ( string .eq. 'Y11') yo(11) = Value
             if ( string .eq. 'Y12') yo(12) = Value
             if ( string .eq. 'Y13') yo(13) = Value
             if ( string .eq. 'Y14') yo(14) = Value
             if ( string .eq. 'Y15') yo(15) = Value
             if ( string .eq. 'Y16') yo(16) = Value
             if ( string .eq. 'Y17') yo(17) = Value
             if ( string .eq. 'Y18') yo(18) = Value
             if ( string .eq. 'Y19') yo(19) = Value
             if ( string .eq. 'Y20') yo(20) = Value

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


