      double complex function A0phiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
C----Expresssion of Eq. (A.4) of arXiv:0906008
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      A0phiqarbmppm=-za(j1,j4)**2/(za(j1,j2)*za(j3,j4))
      return
      end


      double complex function A0phiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
C----Expresssion of Eq. (A.4) of arXiv:0906008
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiqarbmppm
      A0phiqarbmpmp=-A0phiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

      double complex function A0phidqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
C----Expresssion of Eq. (A.4) of arXiv:0906008
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      A0phidqarbmppm=-zb(j2,j3)**2/(zb(j1,j2)*zb(j3,j4))
      return
      end

      double complex function A0phidqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
C----Expresssion of Eq. (A.4) of arXiv:0906008
      include 'constants.f'
      include 'zprods_decl.f'
      double complex A0phidqarbmppm
      integer j1,j2,j3,j4
      A0phidqarbmpmp=-A0phidqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end
