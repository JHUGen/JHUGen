      double complex function Alcphiqarbmpmp(j1,j2,j3,j4,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'deltar.f'
C     arXIv:09060008v1, Eq.4.9
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      double precision s12,s13,s41,s23,s24,s34
      double precision s123,s234,s341,s412,mhsq
      double complex Vlc,Lsm1DS,Lsm1_2me,L0,L1,lnrat,
     & A0phiqarbmpmp,A0phidqarbmpmp
      integer j1,j2,j3,j4

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
      Vlc = -2d0*epinv**2
     &  -epinv*lnrat(musq,-s23)-epinv*lnrat(musq,-s41)
     &  -0.5d0*lnrat(musq,-s23)**2-0.5d0*lnrat(musq,-s41)**2
     & +13d0/6d0*(2d0*epinv+lnrat(musq,-s12)+lnrat(musq,-s34))
     & - Lsm1_2me(s123,s234,s23,mhsq)
     % - Lsm1_2me(s341,s412,s41,mhsq)
     & - Lsm1DS(s23,s34,s234)
     & - Lsm1DS(s34,s41,s341)
     & - Lsm1DS(s41,s12,s412)
     & - Lsm1DS(s12,s23,s123)
     & + 101d0/9d0 - deltar/3d0

C---arXIv:09060008v1 Eq.(4.13)
      Alcphiqarbmpmp=
     & A0phiqarbmpmp(j1,j2,j3,j4,za,zb)*(Vlc
     & +(1d0 -(za(j1,j4)*za(j2,j3)/(za(j1,j3)*za(j2,j4)))**2)
     & *(Lsm1DS(s23,s34,s234) +Lsm1DS(s41,s12,s412)))
     & - 0.5d0*za(j1,j2)*za(j3,j4)/za(j2,j4)**2
     & *(s41**2*L1(-s412,-s12)/s12**2+s23**2*L1(-s234,-s34)/s34**2
     & - lnrat(-s412,-s12)-lnrat(-s234,-s34))
     & + 2*zb(j2,j4)/za(j2,j4)
     & *(za(j2,j3)/zb(j1,j4)
     &  *(s24*L0(-s412,-s12)/s12+lnrat(-s412,-s12))
     & +za(j1,j4)/zb(j2,j3)
     & *(s24*L0(-s234,-s34)/s34+lnrat(-s234,-s34)))
     & + za(j1,j4)*za(j2,j3)*zb(j2,j4)/za(j2,j4)
     &  *(L0(-s412,-s41)/s41+L0(-s234,-s23)/s23)
     & - 0.5d0/za(j2,j4)**2
     & *(za(j3,j4)*(s24-s41)/zb(j1,j2) 
     & +za(j1,j2)*(s24-s23)/zb(j3,j4))
     & -2d0*A0phidqarbmpmp(j1,j2,j3,j4,za,zb)
 
      return
      end


