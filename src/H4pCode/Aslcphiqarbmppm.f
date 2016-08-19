      double complex function Aslcphiqarbmppm(j1,j2,j3,j4,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'deltar.f'
      double precision s12,s13,s41,s23,s24,s34
      double precision s123,s234,s341,s412,mhsq
      double complex Vslc,Lsm1_2me,L0,L1,lnrat,A0phiqarbmppm
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

C---arXIv:09060008v1 Eq.(4.11)
      Vslc=
     & -epinv**2*2d0
     & -epinv*lnrat(musq,-s12)-epinv*lnrat(musq,-s34)
     & -0.5d0*lnrat(musq,-s12)**2-0.5d0*lnrat(musq,-s34)**2
     & -3d0/2d0*(2d0*epinv+lnrat(musq,-s12)+lnrat(musq,-s34))
     & -Lsm1_2me(s412,s123,s12,mhsq)
     & -Lsm1_2me(s234,s341,s34,mhsq)
     & -7d0-deltar

C---arXIv:09060008v1 Eq.(4.10)
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      Aslcphiqarbmppm=A0phiqarbmppm(j1,j2,j3,j4,za,zb)*Vslc
     & +0.5d0*za(j1,j2)*za(j3,j4)*zb(j2,j3)**2
     & *(L1(-s123,-s12)/s12**2 + L1(-s234,-s34)/s34**2)
     & - za(j1,j4)*zb(j2,j3)*(L0(-s123,-s12)/s12+L0(-s234,-s34)/s34)
      return
      end

