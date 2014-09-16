C-----------------------------------------------------------------------
C THIS COMMON BLOCK CONTAINS LOCAL ARRAY ADDRESSES USED BY ASA OF THE
C CONVERTED PACE VERSION FOR STAGGERED GRID
C-----------------------------------------------------------------------


C ADDRESS GROUP 1

      integer ju, jv, jw, jp, jsc1, jsc1var, jsc1dis, jsc1varm, jsc1dism
      integer js, jres, jgam, jrho, jap, jaw, jae, jas, jan, jal, jar
      integer jgamsc1, jgamlam
      integer jqi, jqj, jqk
      integer jx, jy, jz
      integer ju_int, jv_int, jw_int, jhn_u, jhn_v, jhn_w, jhn_sc1
      integer ju_m, jv_m, jw_m, jf_m, jT_m, jp_m, jrho_m
      integer jruu, jrvv, jrww, jruv, jrvw, jrwu, jrff, jrrh
      integer jsuu, jsvv, jsww, jsuv, jsvw, jswu
      integer jsuf, jsvf, jswf
      integer jtw_sfi, jtw_nfi, jtw_lfi, jtw_rfi,
     &        jtw_wfj, jtw_efj, jtw_lfj, jtw_rfj,
     &        jtw_wfk, jtw_efk, jtw_sfk, jtw_nfk
      integer jrhonp1, jrhonm1, jDrhonp1Dt
      integer jTemperature, jYCH4, jYO2, jYH2O, jYCO2, jYN2
      integer jZlib, jZvarlib,jrholib,jTlib
      integer jYCH4lib, jYO2lib, jYH2Olib, jYCO2lib, jYN2lib
      integer jsrhot, jsrhoux
      integer js_update
      integer jYH2Om, jYCO2m

      common /address1/
     &  ju, jv, jw, jp, jsc1, jsc1var, jsc1dis, jsc1varm, jsc1dism,
     &  js, jres, jgam, jrho, jap, jaw, jae, jas, jan, jal, jar,
     &  jgamsc1, jgamlam,
     &  jqi, jqj, jqk,
     &  jx, jy, jz,
     &  ju_int, jv_int, jw_int, jhn_u, jhn_v, jhn_w, jhn_sc1,
     &  ju_m, jv_m, jw_m, jf_m, jT_m,jp_m, jrho_m,
     &  jruu, jrvv, jrww, jruv, jrvw, jrwu, jrff,jrrh,
     &  jsuu, jsvv, jsww, jsuv, jsvw, jswu,
     &  jsuf, jsvf, jswf,
     &  jtw_sfi, jtw_nfi, jtw_lfi, jtw_rfi,
     &  jtw_wfj, jtw_efj, jtw_lfj, jtw_rfj,
     &  jtw_wfk, jtw_efk, jtw_sfk, jtw_nfk,
     &  jrhonp1, jrhonm1, jDrhonp1Dt,
     &  jTemperature, jYCH4, jYO2, jYH2O, jYCO2, jYN2,
     &  jZlib, jZvarlib,jrholib,jTlib,
     &  jYCH4lib, jYO2lib, jYH2Olib, jYCO2lib, jYN2lib,jsrhot,jsrhoux,
     &  js_update, jYH2Om, jYCO2m

C ADDRESS GROUP 2

      integer jdip, jdiu, jdiv, jdix, jrdip, jrdiv, jrdix
      integer jdjp, jdju, jdjv, jdjx, jrdjp, jrdju, jrdjx
      integer jhijp, jhjip, jhkip, jhkjp
      integer jhiju, jhjiu, jhkiu, jhkju
      integer jhijv, jhjiv, jhkiv, jhkjv
      integer jrp, jru, jrv
      integer jdidjp, jrdidjp, jdidju, jrdidju, jdjdiv, jrdjdiv
      integer jdkp, jdkw

c35
      common /address2/
     &  jdip, jdiu, jdiv, jdix, jrdip, jrdiv, jrdix,
     &  jdjp, jdju, jdjv, jdjx, jrdjp, jrdju, jrdjx,
     &  jhijp, jhjip, jhkip, jhkjp,
     &  jhiju, jhjiu, jhkiu, jhkju,
     &  jhijv, jhjiv, jhkiv, jhkjv,
     &  jrp, jru, jrv,
     &  jdidjp, jrdidjp, jdidju, jrdidju, jdjdiv, jrdjdiv,
     &  jdkp, jdkw

C ADDRESS GROUP 3

      integer jgi, jgj, jgk

c3
      common /address3/
     &  jgi, jgj, jgk

C ADDRESS GROUP 4

      integer jwrk
      integer jw1a, jw1b, jw1c, jw1d
      integer jw2a, jw2b, jw2c, jw2d, jw2e, jw2f, jw2g
      integer jw3a, jw3b, jw3c, jw3d, jw3e, jw3f

      integer iwk1a, iwk1b
      integer iwk2a, iwk2b, iwk2c, iwk2d, iwk2e

c20
      common /address4/
     &  jwrk,
     &  jw1a, jw1b, jw1c, jw1d,
     &  jw2a, jw2b, jw2c, jw2d, jw2e, jw2f, jw2g,
     &  jw3a, jw3b, jw3c, jw3d, jw3e, jw3f,
     &  iwk1a, iwk1b, iwk2a, iwk2b, iwk2c, iwk2d, iwk2e

C ADDRESS GROUP 5
C - SPECIAL

      integer jin_prof
      common /address5_1/
     &  jin_prof

      integer jinbcg
      common /address5_2/
     &  jinbcg
