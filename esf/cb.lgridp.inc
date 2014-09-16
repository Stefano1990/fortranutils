C-----------------------------------------------------------------------
C COMMON BLOCK FOR LOCAL GRID PARAMETERS
C-----------------------------------------------------------------------

C - ORGANIZED LIKE THIS FOR EASE OF USE
C    - SHOULD REALLY COME THROUGH ARGUMENT LIST LIKE ARRAYS
C    - ESSENTIAL THESE ARE CHANGED PRIOR TO CALLING "control_block"

C-----------------------------------------------------------------------
C LOCAL GRID DATA (INCLUDED HERE FOR EASE OF USE)
C-----------------------------------------------------------------------

      integer   ni, nj, nk
      integer   nip1, njp1, nkp1
      integer   io, jo, ko
      integer   npts
      
C TEMPORARY PACE ONES - REMOVE LATER

      integer   IGOX, IGOY, JGOX, JGOY

C COMMON
      
      common /cblgridp/ ni, nj, nk,
     &                 nip1, njp1, nkp1,
     &                 io, jo, ko,
     &                 npts,
     &                 igox, igoy, jgox, jgoy

