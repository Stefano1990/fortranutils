C-----------------------------------------------------------------------
C DIMENSIONS FOR THE LES PROGRAM
C-----------------------------------------------------------------------

      integer   maxra
      integer   maxia
      integer   maxca
      integer   maxvar
      integer   maxsv
      integer   maxsbc
      integer   maxpe
      integer   maxmgl
      integer   maxbl
      integer   maxbl2

C     maxra is the large chunk of memory allocated
C     to hold most of the large floating point arrays.
C     You may want to change this (the code will give
C     an estimate of the percentage of the allocated
C     memory that is actually used).
#if defined SGI || defined DEC
c      parameter (maxra    =100000000)  ! suitable for 1 proc case
       parameter (maxra    =53000000)  ! suitable for 1 proc case
c      parameter (maxra    =10000000)  ! suitable for 1 proc case
c      parameter (maxra    = 50000000)  ! suitable for 2 proc case
c      parameter (maxra    = 30000000)  ! suitable for 4 proc case
c      parameter (maxra    = 15000000)  ! suitable for 8 proc case
#elif defined T3E
      parameter (maxra    = 51000000)
c      parameter (maxra    = 50000000)
#endif
      parameter (maxia    =  200000)
      parameter (maxca    =     200)
      parameter (maxvar   =      57)
      parameter (maxsv    =      20)
      parameter (maxsbc   =      20)
      parameter (maxpe    =       1)
      parameter (maxmgl   =       3)
C     Maximum number of blocks in total. It is 
C     important that this be correct.
      parameter (maxbl    =     96) 
C     Maximum number of blocks per rank. To be
C     safe, just make this equal to maxbl; there
C     is very little to be gained (at the moment)
C     from a smaller number (although that may
C     change in the future).
      parameter (maxbl2   =      1)

      integer   log2maxpe

      parameter (log2maxpe = 3)
