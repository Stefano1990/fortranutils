C-----------------------------------------------------------------------
C BOUNDARY CONDITION TYPES FOR FORTRAN PROGRAMS
C-----------------------------------------------------------------------

      integer   BCT_WALL
      integer   BCT_ZERO_GRAD
      integer   BCT_CONVECTIVE
      integer   BCT_FIXED
      integer   BCT_ZG
      integer   BCT_CSP
      integer   BCT_SYMMETRY
      integer   BCT_CENTRE_LINE
      integer   BCT_PERIODIC
      integer   BCT_INTERFACE
      integer   BCT_TRACTION

      parameter (BCT_WALL        =  1)
      parameter (BCT_ZERO_GRAD   =  2)
      parameter (BCT_CONVECTIVE  =  3)
      parameter (BCT_FIXED       =  4)
      parameter (BCT_TRACTION    =  6)
      parameter (BCT_ZG          =  7)
      parameter (BCT_CSP         =  8)
      parameter (BCT_SYMMETRY    =  9)
      parameter (BCT_CENTRE_LINE = 10)
      parameter (BCT_PERIODIC    = 12)
      parameter (BCT_INTERFACE   = 13)
