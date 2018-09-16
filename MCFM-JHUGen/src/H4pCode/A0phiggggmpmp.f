      double complex function A0phiggggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
C----Expresssion of Eq. (2.14) of hep-ph/0804.4149v3
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      A0phiggggmpmp=za(j1,j3)**4
     .             /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))
      return
      end

