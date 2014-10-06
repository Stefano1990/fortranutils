C=======================================================================
C CALCULATE COORDINATES OF COARSE GRID VERTICES
C=======================================================================

      subroutine mg_xyz(coorname,x,work)

      implicit none

C PARAMETERS

#include "dim.inc"
#include "bc_type.inc"
#include "comm.inc"

C INPUT/OUTPUT

      real      x(*), work(*)
      character coorname

C LOCAL

      integer   iblk, ilev
      integer   nip2_c, njp2_c, nkp2_c, nip2_f, njp2_f, nkp2_f
      integer   i_c, j_c, k_c, k_f, ijk_c, ijk_f
      integer   joff_c, joff_f
      integer   imult, jmult, kmult
      integer   poff_c, poff_f

c     integer ierr

C INCLUDE

#include "cb.control.inc"
#include "cb.lgrid.inc"
#include "cb.comm.inc"


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C FILL IN THE OVERLAPING POINTS OF COORDINATE

      do iblk = 1,numbl

        if (coorname.eq.'x' .or. coorname.eq.'y') then
          poff_f = mg_off_2d_xy(iblk,1)
          nip2_f = mg_nijkp1(1,iblk,1) + 1
          njp2_f = mg_nijkp1(2,iblk,1) + 1
          call xyz_extrap(coorname,x(poff_f+1),nip2_f,njp2_f)
        else
          poff_f = mg_off_1d_z(iblk,1)
          nkp2_f = mg_nijkp1(3,iblk,1) + 1
          call xyz_extrap(coorname,x(poff_f+1),nkp2_f,1)
        end if

      end do

      call swap_xyz(coorname,1,x,work)

C CONSTRUCT COARSE GRIDS

      do 10 ilev = 2,nummgl

        do 15 iblk = 1,numbl

C-----------------------------------------------------------------------
C x, y COORDINATES
C-----------------------------------------------------------------------

          if (coorname.eq.'x' .or. coorname.eq.'y') then

C INTERIOR POINTS

            poff_c = mg_off_2d_xy(iblk,ilev  )
            poff_f = mg_off_2d_xy(iblk,ilev-1)

            nip2_c = mg_nijkp1(1,iblk,ilev  ) + 1
            njp2_c = mg_nijkp1(2,iblk,ilev  ) + 1
            nip2_f = mg_nijkp1(1,iblk,ilev-1) + 1
            njp2_f = mg_nijkp1(2,iblk,ilev-1) + 1

            imult = mg_ijk_mult(1,iblk,ilev-1)
            jmult = mg_ijk_mult(2,iblk,ilev-1)

            do j_c = 2,njp2_c-1
              joff_c = poff_c + nip2_c*(j_c-1)
              joff_f = poff_f + nip2_f*(jmult*j_c-3)
              do i_c = 2,nip2_c-1
                ijk_c = joff_c + i_c
                ijk_f = joff_f + imult*i_c-2
                x(ijk_c) = x(ijk_f)
              end do
            end do

          call xyz_extrap(coorname,x(poff_c+1),nip2_c,njp2_c)

C-----------------------------------------------------------------------
C z COORDINATES
C-----------------------------------------------------------------------

          else

C INTERIOR POINTS

            poff_c = mg_off_1d_z(iblk,ilev  )
            poff_f = mg_off_1d_z(iblk,ilev-1)

            nkp2_c = mg_nijkp1(3,iblk,ilev  ) + 1
            nkp2_f = mg_nijkp1(3,iblk,ilev-1) + 1

            kmult = mg_ijk_mult(3,iblk,ilev-1)

            do k_c = 2,nkp2_c-1
              k_f = poff_f + kmult*k_c-2
              x(poff_c+k_c) = x(k_f)
            end do

            call xyz_extrap(coorname,x(poff_c+1),nkp2_c,1)

          end if


 15     continue

C-----------------------------------------------------------------------
C SWAP INTERFACE VALUES
C-----------------------------------------------------------------------

        call swap_xyz(coorname,ilev,x,work)


 10   continue

      return
      end
