      double complex function A0phiggggmmmm(j1,j2,j3,j4,za,zb)
      implicit none
C----Expresssion of Eq. (7) of hep-ph/0607139v2
c--- (or, equivalently, Eq.(2.17) of 0704.3914v3)
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double precision mhsq
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      A0phiggggmmmm=mhsq**2/(zb(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j4,j1))
      return
      end

