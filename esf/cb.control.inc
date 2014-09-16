C-----------------------------------------------------------------------
C GLOBAL CONTROL DATA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C CONTROL
C-----------------------------------------------------------------------

      integer           ncycle(maxvar)

      integer           istep, iststep, numsteps, ststep
      integer           ivar, numvar, numsolve
      integer           psolver

      integer           rotate

      common / cbcontrol1 /
     &                  ncycle, istep, iststep, numsteps, ststep,
     &                  ivar, numvar, numsolve,
     &                  psolver,
     &                  rotate

      real              relax(maxvar)
      real              residm(maxvar)
      real              pref_xyz(3)
      real              thetime ,dtime, phytime

      real              maxDFS

      real              t_relax, alpha

      common / cbcontrol2 /
     &                  relax, residm, pref_xyz, thetime, dtime,
     &                  phytime,
     &                  maxDFS,
     &                  t_relax, alpha

      logical           gotu ,gotv ,gotw, gotp, gotsc1, debugging

      common / cbcontrol3 /
     &                  gotu ,gotv ,gotw, gotp, gotsc1, debugging

      character*32      namvar(maxvar)
      integer           pntvar(maxvar)

      common / cbcontrol4 /
     &                  namvar, pntvar

      integer           bl_ref_block, bl_ref_ijk(3)
      common / cbcontrol5 /
     &                  bl_ref_block, bl_ref_ijk

C-----------------------------------------------------------------------
C INPUT/OUTPUT DEVICE NUMBER
C-----------------------------------------------------------------------

      integer           ndn, nblset
      common / cbcontrol6 /ndn, nblset

C-----------------------------------------------------------------------
C DIAGNOSTIC PRINT
C-----------------------------------------------------------------------

      integer           idiagv(maxbl,maxvar)
      integer           idiag(maxbl)
      integer           idskip(maxbl)
      integer           idchunk(6,maxbl)

      common / cbdiag   /
     &                  idiagv, idiag, idskip, idchunk

C-----------------------------------------------------------------------
C FILE
C-----------------------------------------------------------------------

      integer           rformat, wformat
      integer           ifile(maxvar)
      integer           ftwall, ifskip

      common / cbfile   /
     i                  rformat, wformat, ifile, ftwall, ifskip

C-----------------------------------------------------------------------
C PRINT
C-----------------------------------------------------------------------

      integer           iprint(maxbl,maxvar)
      integer           ipskip(maxbl)
      integer           ipchunk(6,maxbl)

      common / cbdiag   /
     &                  iprint, ipskip, ipchunk

C-----------------------------------------------------------------------
C PROPERTIES
C-----------------------------------------------------------------------

C THIS WANTS CHANGING TO A NON-DIMENSIONAL FORM!!!!

      real              pr(maxvar)

      real              gamma ,rhohi, rholo ,visc, rhoini

      common / cbprop   /
     &                  pr,
     &                  gamma ,rhohi, rholo ,visc, rhoini

C-----------------------------------------------------------------------
C TURBULENCE MODEL
C-----------------------------------------------------------------------

      real              prt(maxvar)

      real              ca    ,ck, aplus, csles


      common / cbturb1 /
     &                  prt,
     &                  ca, ck, aplus, csles

      integer           iturb
      integer           iwf

      common / cbturb2 /
     &                  iturb ,iwf

C-----------------------------------------------------------------------
C GLOBAL MULTIGRID DATA
C-----------------------------------------------------------------------

      integer   meth_rest, meth_prol
      integer   meth_smooth, meth_csolve
      integer   nummgl
      integer   meth_cycle
      integer   nummgits
      integer   numcits, numpre, numpost
      integer   meth_coarse

c---- SMOOTHING CONTROL PARAMETERS FOR FANCY M.G.

CC      real  eps(maxmgl), delta, eta

      common /cbmgdata1/
     &       meth_rest, meth_prol,
     &       meth_smooth, meth_csolve,
     &       nummgl,
     &       meth_cycle,
     &       nummgits,
     &       numcits, numpre, numpost,
     &       meth_coarse

      real      tolc
      real      tolf
      real      relax_pmg
      real      relax_pmgc
      real      rnorm(maxmgl), fnorm(maxmgl), enorm(maxmgl),
     &          vnorm(maxmgl)

      common /cbmgdata2/
     &       tolc,
     &       tolf,
     &       relax_pmg,
     &       relax_pmgc,
     &       rnorm, fnorm, enorm, vnorm

C-----------------------------------------------------------------------
C LEGNTHS OF WORKSPACES
C-----------------------------------------------------------------------

      integer lrwk, lrwk2d, liwk
      common  /cbwrksp/
     &        lrwk, lrwk2d, liwk

C-----------------------------------------------------------------------
C INLET BC GATHERING
C-----------------------------------------------------------------------

      integer inbcglength
      integer numinbcg, inbcgblk, inbcgplane, inbcgfileskip,
     &        inbcgfileoff
      integer inbcgvar(maxvar)

      common /cb_inbcg/
     &       inbcglength,
     &       numinbcg, inbcgblk, inbcgplane, inbcgfileskip,
     &       inbcgfileoff,
     &       inbcgvar

C-----------------------------------------------------------------------
C INLET BC SPECIFICATION
C-----------------------------------------------------------------------

      integer inbcsblk, inbcsplane

      common /cb_inbcs1/
     &       inbcsblk, inbcsplane

      real   inbcsmflux, inbcsperturb(3)
      common /cb_intcs2/
     &       inbcsmflux, inbcsperturb

C-----------------------------------------------------------------------
C STATISTICS RELATED
C-----------------------------------------------------------------------

      integer addback

      common  /cb_statistics/
     &        addback
C----------------------------------------------------------------------
C WAM additions
C----------------------------------------------------------------------
      integer  Stat_Freq, Fiel_Freq, Dump_Freq, Spec_Freq
      logical  New_Run, specout, writewall
      integer  last,step_no,jsplit,no_sam,statstart
      integer  refpresblock, inlet1, inlet2,inlet3, inlet4
      integer  fieldstart 
 
      common /cb_wamadd/
     &         Stat_Freq, Fiel_Freq, Dump_Freq, Spec_Freq,
     &         New_Run,last,step_no,jsplit,no_sam,statstart,
     &         inlet1,inlet2, specout,writewall,inlet3, inlet4,
     &         fieldstart

      real     yz(1:20),xo(1:20),yplane,yo(1:20)
      integer  xoilow(1:20,2),yojlow(1:20,2)

      common /cb_wamadd2/
     &        yz,xo,yo,yplane,refpresblock,
     &        xoilow,yojlow

      logical LES
      integer None, Smagorinsky, WALE, SIGMA,Model
         PARAMETER (None=0, Smagorinsky=1, WALE=2,SIGMA=3)
         
      common /cb_wamadd3/
     & Model
