C-----------------------------------------------------------------------
C COMMON BLOCKS FOR LOCAL BLOCK/GRID DATA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C BLOCKS - LOCAL
C-----------------------------------------------------------------------

      integer      numbl

C LINKS TO OTHER BLOCKS
      
      integer      bl2bl(6,maxbl)
      integer      bl2side(6,maxbl)
      integer      bl2pe(maxbl)
c     integer      ibl2iblt(maxbl)

C BOUNDARY CONDITIONS

      integer      bl_bc(6,maxbl)
      real         bl_split(6,maxbl)

C GRID PARAMETERS

      integer      bl_nijkp1(3,maxbl)
      integer      bl_nijkp1_master(3,maxbl)

C BLOCK OFF SET FOR SINGLE GRID

      integer      sg_off_1d_i(maxbl)
      integer      sg_off_1d_j(maxbl)
      integer      sg_off_1d_k(maxbl)
      integer      sg_off_2d_ij(maxbl)
      integer      sg_off_2d_jk(maxbl)
      integer      sg_off_2d_ki(maxbl)
      integer      sg_off(maxbl)

C MULTIGRID DATA

C CHECK DIMS ARE REALLY NEEDED - LOCAL? ?????????????????

      integer      mg_nijkp1(3,maxbl,maxmgl)
      integer      mg_dim1(maxbl)
      integer      mg_dim2(maxbl)
      integer      mg_dim3(maxbl)
      integer      mg_off_1d_i(maxbl,maxmgl)
      integer      mg_off_1d_j(maxbl,maxmgl)
      integer      mg_off_1d_k(maxbl,maxmgl)
      integer      mg_off_2d_ij(maxbl,maxmgl)
      integer      mg_off_2d_jk(maxbl,maxmgl)
      integer      mg_off_2d_ki(maxbl,maxmgl)
      integer      mg_off(maxbl,maxmgl)
      integer      mg_ijk_mult(3,maxbl,maxmgl)

      integer      mg_off_1d_z(maxbl,maxmgl)
      integer      mg_off_2d_xy(maxbl,maxmgl)
      
      common /cblgrid1_1/
     &       numbl, bl2bl, bl2side, bl2pe, bl_bc, bl_nijkp1,
     &       bl_nijkp1_master

      common /cblgrid1_2/
     &       bl_split

      common /cblgrid2/
     &       sg_off_1d_i, sg_off_1d_j, sg_off_1d_k, sg_off_2d_ij,
     &       sg_off_2d_jk, sg_off_2d_ki, sg_off,
     &       mg_nijkp1, mg_dim1, mg_dim2, mg_dim3,
     &       mg_off_1d_i, mg_off_1d_j, mg_off_1d_k,
     &       mg_off_2d_ij, mg_off_2d_jk, mg_off_2d_ki,
     &       mg_off, mg_ijk_mult,
     &       mg_off_1d_z, mg_off_2d_xy

      integer      maxijk1, maxijk2, maxijk3
      integer      maxijk1_mg, maxijk2_mg, maxijk3_mg

      common /cblgrid3/
     &       maxijk1, maxijk2, maxijk3,
     &       maxijk1_mg, maxijk2_mg, maxijk3_mg

C-----------------------------------------------------------------------
C LENGTH PARAMETERS RELATED TO GRID SIZE
C-----------------------------------------------------------------------

      integer bl_tlen_mg_i, bl_tlen_mg_j, bl_tlen_mg_k, bl_tlen_mg_ij,
     &        bl_tlen_mg_ijk
      integer bl_tlen_sg_ij, bl_tlen_sg_jk, bl_tlen_sg_ki,
     &        bl_tlen_sg_ijk
      integer bl_tlen_mg_z, bl_tlen_mg_xy
      common  /cblgrid4/
     &        bl_tlen_mg_i, bl_tlen_mg_j, bl_tlen_mg_k, bl_tlen_mg_ij,
     &        bl_tlen_mg_ijk,
     &        bl_tlen_sg_ij, bl_tlen_sg_jk, bl_tlen_sg_ki,
     &        bl_tlen_sg_ijk,
     &        bl_tlen_mg_z, bl_tlen_mg_xy

C-----------------------------------------------------------------------
C STORAGE FOR MASS FLUXES FOR CONVECTIVE BOUNDARIES
C-----------------------------------------------------------------------

      real    mfluxbl_out(6,maxbl), mfluxbl_out2(6,maxbl),
     &        area_exitbl(6,maxbl),area_rho_exitbl(6,maxbl),
     &        mfluxbl_out_ttl(6), area_rho_exitbl_ttl(6)
      common  /cblgrid5/
     &        mfluxbl_out, mfluxbl_out2, area_exitbl,area_rho_exitbl,
     &        mfluxbl_out_ttl, area_rho_exitbl_ttl
