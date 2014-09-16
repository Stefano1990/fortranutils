      integer imax, jmax, kmax
      parameter (imax=66, jmax=258, kmax=303)
      real    dfdxw(maxbl,jmax,kmax), dfdxe(maxbl,jmax,kmax),
     &        dfdys(maxbl,imax,kmax), dfdyn(maxbl,imax,kmax),
     &        dfdzl(maxbl,imax,jmax), dfdzr(maxbl,imax,jmax)
      common /grad/  dfdxw, dfdxe, dfdys, dfdyn, dfdzl, dfdzr
