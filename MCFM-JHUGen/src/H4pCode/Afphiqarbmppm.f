      double complex function Afphiqarbmppm(j1,j2,j3,j4,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      double precision s12,s34
      double complex lnrat,A0phiqarbmppm
      integer j1,j2,j3,j4
      s12=s(j1,j2)
      s34=s(j3,j4)
C---arXIv:09060008v1 Eq.(4.12)
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      Afphiqarbmppm=A0phiqarbmppm(j1,j2,j3,j4,za,zb)*
     & (-2d0/3d0*(2d0*epinv+lnrat(musq,-s12)+lnrat(musq,-s34))
     & - 20D0/9D0)

      return
      end
