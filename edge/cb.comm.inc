C-----------------------------------------------------------------------
C DATA FOR COMMUNICATIONS
C-----------------------------------------------------------------------
      
      integer   pe
      integer   numpe
      integer   numbl_total

      integer   rlext
      integer   workwin

      integer   idwrcount
      integer   iccount
      integer   ilrwk
      
      integer   inwoff(6,maxbl,maxmgl)
      integer   inikoff(6,6,maxbl,maxmgl)
      integer   iniknum(2,6,maxbl,maxmgl)
      integer   inoffmax(maxmgl)
      integer   bl_rb(maxbl)

C     These aren't really necessary. But I still think they ought
C     to be there for convenience.
      integer   iblkg2l(maxbl)
      integer   iblkl2g(maxbl)

      logical   master

      common / cbcomm / pe, numpe, numbl_total,
     &                  inwoff, inoffmax, inikoff, iniknum,
     &                  bl_rb,
     &                  iblkg2l, iblkl2g, rlext, workwin,
     &                  idwrcount, iccount, ilrwk, master
      

