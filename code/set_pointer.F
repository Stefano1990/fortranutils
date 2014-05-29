C=======================================================================
C SPLIT THE LARGE ARRAYS AND SET UP THE POINTERS OF ALL ARRAY QUANTITIES
C=======================================================================

      subroutine set_pointer()

      implicit none

C PARAMETERS

#include "dim.inc"

C LOCAL

      integer ii, ir
c     integer bl_tlen_mg_surf

C INCLUDE

#include "cb.control.inc"
#include "cb.address.inc"
#include "cb.lgrid.inc"

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C REAL ARRAYS
C-----------------------------------------------------------------------

C WORKSPACE


C - 1-D WORKSPACE

      ir = 1

      jw1a = ir
      jw1b = ir + maxijk1
      jw1c = ir + maxijk1*2
      jw1d = ir + maxijk1*3
      ir   = ir + maxijk1*4

C - 2-D WORKSPACE

      jw2a = ir
      jw2b = ir + maxijk2
      jw2c = ir + maxijk2*2
      jw2d = ir + maxijk2*3
      jw2e = ir + maxijk2*4
      jw2f = ir + maxijk2*5
      jw2g = ir + maxijk2*6
      ir   = ir + maxijk2
      ir   = ir + maxijk2*MAX(MAXBL2+1,7)
      lrwk2d= ir - jw1a +1


C - 3-D WORKSPACE

      jw3a = ir
      jw3b = ir + maxijk3
      jw3c = ir + maxijk3*2
      jw3d = ir + maxijk3*3
      jw3e = ir + maxijk3*4
      jw3f = ir + maxijk3*5

c      ir = ir + maxijk3*6
      ir = ir + maxijk3

      ir = ir + MAXIJK3*max(maxbl2+1,6)

      lrwk = ir - jw1a + 1

C - COORDINATES AND METRIC QUANTITIES

      bl_tlen_mg_i  =
     &  mg_off_1d_i(numbl,nummgl) + mg_nijkp1(1,numbl,nummgl)
      bl_tlen_mg_j  =
     &  mg_off_1d_j(numbl,nummgl) + mg_nijkp1(2,numbl,nummgl)
      bl_tlen_mg_k  =
     &  mg_off_1d_k(numbl,nummgl) + mg_nijkp1(3,numbl,nummgl)
      bl_tlen_mg_ij = mg_off_2d_ij(numbl,nummgl) +
     &  mg_nijkp1(1,numbl,nummgl)*mg_nijkp1(2,numbl,nummgl)

      bl_tlen_mg_z  =
     &  mg_off_1d_z(numbl,nummgl) + mg_nijkp1(3,numbl,nummgl)+1
      bl_tlen_mg_xy = mg_off_2d_xy(numbl,nummgl) +
     &  (mg_nijkp1(1,numbl,nummgl)+1)*(mg_nijkp1(2,numbl,nummgl)+1)

C -- x, y COORDINATES

      jx = ir
      jy = ir + bl_tlen_mg_xy

      ir = ir + bl_tlen_mg_xy*2

C -- z COORDINATE

      jz = ir
      ir = ir + bl_tlen_mg_z

C -- METRIC QUANTTITIES

      jdkp  = ir
      jdkw  = ir + bl_tlen_mg_k

      ir = ir + bl_tlen_mg_k*2

      jdip    = ir
      jdiu    = ir + bl_tlen_mg_ij
      jdiv    = ir + bl_tlen_mg_ij*2
      jdix    = ir + bl_tlen_mg_ij*3
      jrdip   = ir + bl_tlen_mg_ij*4
      jrdiv   = ir + bl_tlen_mg_ij*5
      jrdix   = ir + bl_tlen_mg_ij*6
      jdjp    = ir + bl_tlen_mg_ij*7
      jdju    = ir + bl_tlen_mg_ij*8
      jdjv    = ir + bl_tlen_mg_ij*9
      jdjx    = ir + bl_tlen_mg_ij*10
      jrdjp   = ir + bl_tlen_mg_ij*11
      jrdju   = ir + bl_tlen_mg_ij*12
      jrdjx   = ir + bl_tlen_mg_ij*13
      jhijp   = ir + bl_tlen_mg_ij*14
      jhjip   = ir + bl_tlen_mg_ij*15
      jhkip   = ir + bl_tlen_mg_ij*16
      jhkjp   = ir + bl_tlen_mg_ij*17
      jhiju   = ir + bl_tlen_mg_ij*18
      jhjiu   = ir + bl_tlen_mg_ij*19
      jhkiu   = ir + bl_tlen_mg_ij*20
      jhkju   = ir + bl_tlen_mg_ij*21
      jhijv   = ir + bl_tlen_mg_ij*22
      jhjiv   = ir + bl_tlen_mg_ij*23
      jhkiv   = ir + bl_tlen_mg_ij*24
      jhkjv   = ir + bl_tlen_mg_ij*25
      jrp     = ir + bl_tlen_mg_ij*26
      jru     = ir + bl_tlen_mg_ij*27
      jrv     = ir + bl_tlen_mg_ij*28
      jdidjp  = ir + bl_tlen_mg_ij*29
      jrdidjp = ir + bl_tlen_mg_ij*30
      jdidju  = ir + bl_tlen_mg_ij*31
      jrdidju = ir + bl_tlen_mg_ij*32
      jdjdiv  = ir + bl_tlen_mg_ij*33
      jrdjdiv = ir + bl_tlen_mg_ij*34

      ir = ir + bl_tlen_mg_ij*35

C - SINGLE GRID QUANTITIES

      bl_tlen_sg_ijk = sg_off(numbl) +
     &  bl_nijkp1(1,numbl)*bl_nijkp1(2,numbl)*bl_nijkp1(3,numbl)

C -- VARIABLES

      ju = ir
      jv = ir + bl_tlen_sg_ijk
      jw = ir + bl_tlen_sg_ijk*2

      ju_int = ir + bl_tlen_sg_ijk*3
      jv_int = ir + bl_tlen_sg_ijk*4
      jw_int = ir + bl_tlen_sg_ijk*5

      jhn_u = ir + bl_tlen_sg_ijk*6
      jhn_v = ir + bl_tlen_sg_ijk*7
      jhn_w = ir + bl_tlen_sg_ijk*8

      ir = ir + bl_tlen_sg_ijk*9

      jsc1 = ir
      jhn_sc1 = ir + bl_tlen_sg_ijk

      ir = ir + bl_tlen_sg_ijk*2

      jsc1var = ir
      jsc1dis = ir + bl_tlen_sg_ijk
      ir = ir + bl_tlen_sg_ijk*2

C -- GAM

      jgam = ir
      ir = ir + bl_tlen_sg_ijk
      jgamsc1 = ir
      ir = ir + bl_tlen_sg_ijk
      jgamlam = ir
      ir = ir + bl_tlen_sg_ijk

C -- BOUNDARY VALUE STORAGE (TEMPORARILY NOT USED. BUT FOR GENERAL CASE,
C    THEY SHOULD BE USED)

c     bl_tlen_mg_surf = mg_off_mg_surf(numbl,nummgl) +
c    &  2*mg_nijkp1(1,numbl,nummgl)*mg_nijkp1(2,numbl,nummgl) +
c    &  2*mg_nijkp1(2,numbl,nummgl)*mg_nijkp1(3,numbl,nummgl) +
c    &  2*mg_nijkp1(3,numbl,nummgl)*mg_nijkp1(1,numbl,nummgl)

c     ju_bv   = ir
c     jv_bv   = ir + bl_tlen_mg_surf
c     jw_bv   = ir + bl_tlen_mg_surf*2
c     jp_bv   = ir + bl_tlen_mg_surf*3
c     jgam_bv = ir + bl_tlen_mg_surf*4

c     ir = ir + bl_tlen_mg_surf*5

C -- STATISTICS

      bl_tlen_sg_ij = sg_off_2d_ij(numbl) +
     &                bl_nijkp1(1,numbl)*bl_nijkp1(2,numbl)

      ju_m = ir
      jv_m = ir + bl_tlen_sg_ij
      jw_m = ir + bl_tlen_sg_ij*2

      ir = ir + bl_tlen_sg_ij*3

      jf_m = ir
      ir = ir + bl_tlen_sg_ij

      jsc1varm = ir
      ir = ir + bl_tlen_sg_ij
      jsc1dism = ir

      jruu = ir + bl_tlen_sg_ij
      jrvv = ir + bl_tlen_sg_ij*2
      jrww = ir + bl_tlen_sg_ij*3
      jruv = ir + bl_tlen_sg_ij*4
      jrvw = ir + bl_tlen_sg_ij*5
      jrwu = ir + bl_tlen_sg_ij*6
      jrff = ir + bl_tlen_sg_ij*7
      jrrh = ir + bl_tlen_sg_ij*8
      jp_m = ir + bl_tlen_sg_ij*9
      jrho_m = ir + bl_tlen_sg_ij*10
      
      jsuu = ir + bl_tlen_sg_ij*11
      jsvv = ir + bl_tlen_sg_ij*12
      jsww = ir + bl_tlen_sg_ij*13
      jsuv = ir + bl_tlen_sg_ij*14
      jsvw = ir + bl_tlen_sg_ij*15
      jswu = ir + bl_tlen_sg_ij*16

      ir = ir + bl_tlen_sg_ij*17

C -- WALL SHEAR STRESS

      bl_tlen_sg_ij = sg_off_2d_ij(numbl) +
     &                bl_nijkp1(1,numbl)*bl_nijkp1(2,numbl)
      bl_tlen_sg_jk = sg_off_2d_jk(numbl) +
     &                bl_nijkp1(2,numbl)*bl_nijkp1(3,numbl)
      bl_tlen_sg_ki = sg_off_2d_ki(numbl) +
     &                bl_nijkp1(3,numbl)*bl_nijkp1(1,numbl)

      jtw_sfi = ir
      jtw_nfi = jtw_sfi + bl_tlen_sg_ki
      jtw_lfi = jtw_nfi + bl_tlen_sg_ki
      jtw_rfi = jtw_lfi + bl_tlen_sg_ij

      jtw_wfj = jtw_rfi + bl_tlen_sg_ij
      jtw_efj = jtw_wfj + bl_tlen_sg_jk
      jtw_lfj = jtw_efj + bl_tlen_sg_jk
      jtw_rfj = jtw_lfj + bl_tlen_sg_ij

      jtw_wfk = jtw_rfj + bl_tlen_sg_ij
      jtw_efk = jtw_wfk + bl_tlen_sg_jk
      jtw_sfk = jtw_efk + bl_tlen_sg_jk
      jtw_nfk = jtw_sfk + bl_tlen_sg_ki

      ir = jtw_nfk + bl_tlen_sg_ki

C - MULTIGRID QUANTITIES

      bl_tlen_mg_ijk = mg_off(numbl,nummgl) +
     &  mg_nijkp1(1,numbl,nummgl)*mg_nijkp1(2,numbl,nummgl)*
     &  mg_nijkp1(3,numbl,nummgl)

C -- DENSITY

      jrho = ir
      ir = ir + bl_tlen_mg_ijk

      jrhonp1 = ir
      ir = ir + bl_tlen_mg_ijk

      jrhonm1 = ir
      ir = ir + bl_tlen_mg_ijk

      jDrhonp1Dt = ir
      ir = ir + bl_tlen_mg_ijk

      jTemperature = ir
      ir = ir + bl_tlen_mg_ijk

      jT_m = ir
      ir = ir + bl_tlen_mg_ijk


C -- VARIABLES

      jp = ir
      ir = ir + bl_tlen_mg_ijk

C -- MATRIX ENTRIES, SOURCE TERM, AND RESIDUAL

      jap = ir
      jae = ir + bl_tlen_mg_ijk
      jaw = ir + bl_tlen_mg_ijk*2
      jas = ir + bl_tlen_mg_ijk*3
      jan = ir + bl_tlen_mg_ijk*4
      jal = ir + bl_tlen_mg_ijk*5
      jar = ir + bl_tlen_mg_ijk*6

      js = ir + bl_tlen_mg_ijk*7

      jres = ir + bl_tlen_mg_ijk*8

      ir = ir + bl_tlen_mg_ijk*9

      jsrhot = ir
      ir = ir + bl_tlen_mg_ijk

      jsrhoux = ir
      ir = ir + bl_tlen_mg_ijk

      js_update = ir
      ir = ir + bl_tlen_mg_ijk


C -- SMOOTHER RELEVANT ARRAYS

C NOT EXTRA ARRAYS INVOLVED FOR USING POINTER SMOOTHER

C -- ARRAYS STORING THE VARIALBES ON A PLANE WHICH ARE TO BE USED AS
C    BOUNDARY CONDITIONS

      jinbcg = ir
      ir = ir + inbcglength

C -- ARRAYS STORING SPECIFIED VELOCITY PROFILES (ASSUME ONLY FOR
C    u,v,w)

      jin_prof = ir
      ir = ir  + bl_tlen_mg_k*3

C-----------------------------------------------------------------------
C INTEGER ARRAYS
C-----------------------------------------------------------------------

      ii = 1

      iwk1a = ii

      liwk = maxia - iwk1a + 1

C-----------------------------------------------------------------------
C CHECK
C-----------------------------------------------------------------------

      if (ir .gt. maxra) then
        write(*,'(''Real array size requires adjustment: '',
     &    ''Declared: '',i,''; Required: '',i)') maxra, ir
        stop
      else
        write(*,'(''Real array length: '',i20,
     &    '', '',f5.2,''% used'')') maxra,float(ir)/float(maxra)*100.
      end if
      if (ir .gt. maxra) then
        write(*,'(''Integer array size requires adjustment: '',
     &    ''Length: '',i,''; Required: '',i)') maxia, ii
        stop
      else
        write(*,'(''Integer array length: '',i20,
     &    '', '',f5.2,''% used'')') maxia,float(ii)/float(maxia)*100.
      end if

      return
      end
