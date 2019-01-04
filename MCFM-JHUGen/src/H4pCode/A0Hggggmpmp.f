      double complex function A0Hggggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiggggmpmp
      A0Hggggmpmp=A0phiggggmpmp(j1,j2,j3,j4,za,zb)
     .           +A0phiggggmpmp(j4,j1,j2,j3,zb,za)
      return
      end
