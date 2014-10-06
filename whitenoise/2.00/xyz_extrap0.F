C=======================================================================
C EXTRAPOLATE x, y, z COORDINATES IN ORDER TO ALLOW VIRTUAL CELLS TO BE
C COVERED
C=======================================================================

C THIS IS NOT ACCURATE AT ALL IF CELL ASPECT RATIO CHANGES TOO MUCH

      subroutine xyz_extrap(coorname,x,nip2,njp2)

      implicit  none

C INPUT/OUTPUT

      integer   nip2, njp2
      real      x(*)
      character coorname

C LOCAL

      integer   joff
      integer   ni, nj, nip1, njp1
      integer   i, j, ij

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      nip1 = nip2 - 1
      ni   = nip2 - 2

      if (coorname.eq.'x' .or. coorname.eq.'y') then

        njp1 = njp2 - 1
        nj   = njp2 - 2

        do j = 1,njp2
          joff = nip2*(j-1)
          ij = 2    + joff
          x(ij-1) = 2.0*x(ij) - x(ij+1)
          ij = nip1 + joff
          x(ij+1) = 2.0*x(ij) - x(ij-1)
        end do

        joff = nip2*nj
        do i = 1,nip2
          ij = i + nip2
          x(ij-nip2) = 2.0*x(ij) - x(ij+nip2)
          ij = i + joff
          x(ij+nip2) = 2.0*x(ij) - x(ij-nip2)
        end do

      else

        x(1   ) = 2.0*x(2   ) - x(3 )
        x(nip2) = 2.0*x(nip1) - x(ni)

      end if

      return
      end


