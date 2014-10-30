C-----------------------------------------------------------------------
C COMMON BLOCKS FOR MASTER BLOCK/GRID DATA
C-----------------------------------------------------------------------

C NOTE: THIS IS ONLY USED BY THE MASTER PROCESSOR TO ASSEMBLE DATA FOR
C       THE LOCAL DATA FOR EACH PROCESSOR

C-----------------------------------------------------------------------
C BLOCKS
C-----------------------------------------------------------------------

      integer      nummbl
      
      integer      mbl_nijkp1(3,maxmbl)
      integer      mbl_2_mbl(6,maxmbl)
      integer      mbl_2_side(6,maxmbl)
      integer      mbl_pe(maxmbl)

C BOUNDARY CONDITIONS

      integer      mbl_bc(6,maxmbl)
      real         mbl_split(6,maxmbl)

C PRESSURE REFERERENCE POINT

      integer      mbl_ref_block
      integer      mbl_ref_ijk

C WORK ARRAYS - MUST BE SAME SIZE ON ALL PROCESSORS FOR "shmem"

      integer      mbl_max_nijkp1(3)
      integer      mbl_maxijk1
      integer      mbl_maxijk2
      integer      mbl_maxijk3

      common /cbcomm1_1/
     &       nummbl, mbl_nijkp1, mbl_2_mbl, mbl_2_side, mbl_pe,
     &       mbl_bc, mbl_ref_block, bml_ref_ijk,
     &       mbl_max_nijkp1, mbl_maxijk1, mbl_maxijk2, mbl_maxijk3

      common /cbcomm1_1/
     &        mbl_split
