      double complex function A0Hggggmmpp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiggggmmpp
      A0Hggggmmpp=A0phiggggmmpp(j1,j2,j3,j4,za,zb)
     .           +A0phiggggmmpp(j3,j4,j1,j2,zb,za)
      return
      end
