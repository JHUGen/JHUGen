      function Alcphiqarbmpmp(j1,j2,j3,j4,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Alcphiqarbmpmp
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'deltar.f'
C     arXIv:09060008v1, Eq.4.9
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      real(dp):: s12,s13,s41,s23,s24,s34
      real(dp):: s123,s234,s341,s412,mhsq
      complex(dp):: Vlc,Lsm1DS,Lsm1_2me,L0,L1,lnrat,
     & A0phiqarbmpmp,A0phidqarbmpmp
      integer:: j1,j2,j3,j4

      s12=s(j1,j2)
      s13=s(j1,j3)
      s41=s(j4,j1)
      s23=s(j2,j3)
      s24=s(j2,j4)
      s34=s(j3,j4)
      s123=s12+s13+s23
      s234=s23+s24+s34
      s341=s34+s41+s13
      s412=s41+s12+s24
      mhsq=s12+s13+s41+s23+s24+s34

C---arXIv:09060008v1 Eq.(4.9)
      Vlc = -2._dp*epinv**2
     &  -epinv*lnrat(musq,-s23)-epinv*lnrat(musq,-s41)
     &  -0.5_dp*lnrat(musq,-s23)**2-0.5_dp*lnrat(musq,-s41)**2
     & +13._dp/6._dp*(2._dp*epinv+lnrat(musq,-s12)+lnrat(musq,-s34))
     & - Lsm1_2me(s123,s234,s23,mhsq)
     % - Lsm1_2me(s341,s412,s41,mhsq)
     & - Lsm1DS(s23,s34,s234)
     & - Lsm1DS(s34,s41,s341)
     & - Lsm1DS(s41,s12,s412)
     & - Lsm1DS(s12,s23,s123)
     & + 101._dp/9._dp - deltar/3._dp

C---arXIv:09060008v1 Eq.(4.13)
      Alcphiqarbmpmp=
     & A0phiqarbmpmp(j1,j2,j3,j4,za,zb)*(Vlc
     & +(1._dp -(za(j1,j4)*za(j2,j3)/(za(j1,j3)*za(j2,j4)))**2)
     & *(Lsm1DS(s23,s34,s234) +Lsm1DS(s41,s12,s412)))
     & - 0.5_dp*za(j1,j2)*za(j3,j4)/za(j2,j4)**2
     & *(s41**2*L1(-s412,-s12)/s12**2+s23**2*L1(-s234,-s34)/s34**2
     & - lnrat(-s412,-s12)-lnrat(-s234,-s34))
     & + 2*zb(j2,j4)/za(j2,j4)
     & *(za(j2,j3)/zb(j1,j4)
     &  *(s24*L0(-s412,-s12)/s12+lnrat(-s412,-s12))
     & +za(j1,j4)/zb(j2,j3)
     & *(s24*L0(-s234,-s34)/s34+lnrat(-s234,-s34)))
     & + za(j1,j4)*za(j2,j3)*zb(j2,j4)/za(j2,j4)
     &  *(L0(-s412,-s41)/s41+L0(-s234,-s23)/s23)
     & - 0.5_dp/za(j2,j4)**2
     & *(za(j3,j4)*(s24-s41)/zb(j1,j2) 
     & +za(j1,j2)*(s24-s23)/zb(j3,j4))
     & -2._dp*A0phidqarbmpmp(j1,j2,j3,j4,za,zb)
 
      return
      end


