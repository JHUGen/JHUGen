      function A1ggggmmmmCC(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1ggggmmmmCC

C--- Expresssion of Eq. (19) of Badger and Glover, hep-ph/0607139v2
c--- Note that the corresponding phigggg-dagger amplitude is zero (Eq.(24))
c--- so that this is the full cut-constructible result for the Higgs
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'scprods_com.f'
      include 'zprods_decl.f'
      complex(dp):: A0phiggggmmmm,F31m,F42me,F41m
      integer:: j1,j2,j3,j4,i
      integer,parameter::ii(7)=(/1,2,3,4,1,2,3/)

c--- set up 's-comma' products
      sc(1,1)=zip
      sc(1,2)=s(j1,j2)
      sc(1,3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      sc(1,4)=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,1)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,2)=zip
      sc(2,3)=s(j2,j3)
      sc(2,4)=s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,1)=s(j3,j4)+s(j3,j1)+s(j4,j1)
      sc(3,2)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,3)=zip
      sc(3,4)=s(j3,j4)
      sc(4,1)=s(j4,j1)
      sc(4,2)=s(j4,j1)+s(j4,j2)+s(j1,j2)
      sc(4,3)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(4,4)=zip

      A1ggggmmmmCC=czip
      do i=1,4
c--- NOTE: Arguments of F41m , F42m have been changed to be consistent
c---       with later papers: F41m(s,t;Psq) -> F41m(Psq;s,t)
c---                          F42me(s,t;Psq,Qsq) -> F42me(Psq,Qsq;s,t)
      A1ggggmmmmCC=A1ggggmmmmCC
     &  +F31m(sc(ii(i),ii(i+2)))-F31m(sc(ii(i),ii(i+3)))
     & -0.5_dp*F42me(sc(ii(i),ii(i+3)),sc(ii(i+1),ii(i+2)),
     &              sc(ii(i),ii(i+2)),sc(ii(i+1),ii(i+3)))
     & -0.5_dp*F41m(sc(ii(i),ii(i+2)),sc(ii(i),ii(i+1)),
     &                               sc(ii(i+1),ii(i+2)))
      enddo
      A1ggggmmmmCC=A0phiggggmmmm(j1,j2,j3,j4,za,zb)*A1ggggmmmmCC
      return
      end

      function A1ggggmmmmNCC(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1ggggmmmmNCC

C--- Expresssion of sum of permutations in
c---      Eq. (25) of Badger and Glover, hep-ph/0607139v2
c--- Note that this is the expression for the full Higgs
c--- non cut-constructible part
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      real(dp):: Np
      complex(dp):: A1ggggmmmmsub
c--- Note: overall factor of (fourpi)**2 removed, c.f. Eq. (23)
      Np=2._dp*(1._dp-real(nflav,dp)/xn)
      A1ggggmmmmNCC=Np/6._dp
     & *(+A1ggggmmmmsub(j1,j2,j3,j4,za,zb)
     &   +A1ggggmmmmsub(j2,j3,j4,j1,za,zb)
     &   +A1ggggmmmmsub(j3,j4,j1,j2,za,zb)
     &   +A1ggggmmmmsub(j4,j1,j2,j3,za,zb))
      return
      end


      function A1ggggmmmmsub(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1ggggmmmmsub

C--- Expresssion of Eq. (25) of Badger and Glover, hep-ph/0607139v2
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      real(dp):: s3
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      A1ggggmmmmsub=-s(j1,j3)*zab2(j4,j1,j3,j2)**2
     & /(s3(j1,j2,j3)*zb(j1,j2)**2*zb(j2,j3)**2)
     & +(za(j3,j4)/zb(j1,j2))**2
     & +2._dp*za(j3,j4)*za(j4,j1)/(zb(j1,j2)*zb(j2,j3))
     & +(s(j1,j2)*s(j3,j4)+s3(j1,j2,j3)*s3(j2,j3,j4)-s(j1,j2)**2)
     & /(2._dp*zb(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j4,j1))
c      A1ggggmmmmsub=A1ggggmmmmsub
      return
      end
