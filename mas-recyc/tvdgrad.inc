      integer imax, jmax, kmax
      parameter (imax=38, jmax=68, kmax=260)
      real    dfdxw(maxbl,jmax,kmax), dfdxe(maxbl,jmax,kmax),
     &        dfdys(maxbl,imax,kmax), dfdyn(maxbl,imax,kmax),
     &        dfdzl(maxbl,imax,jmax), dfdzr(maxbl,imax,jmax)
      common /grad/  dfdxw, dfdxe, dfdys, dfdyn, dfdzl, dfdzr
