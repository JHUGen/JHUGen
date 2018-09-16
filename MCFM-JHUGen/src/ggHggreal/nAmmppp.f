      function nAmmppp(j4,j5,j1,j2,j3)
c--- Calculation of the amplitudes using results of S. Badger
c--- as adapted from routines written by F. Caola

      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::nAmmppp
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::j1,j2,j3,j4,j5
      real(dp)::sH1,sH2,sH3,s125,s145,s234,s345,s3,s4,qsq
      complex(dp)::zab_2_pH3_4,zab_5_pH4_3,zab_2_34_5,zab_4_23_1,
     & zab_2_pH_1,zab_2_pH_3,zbb_1_pH_45_3, zbb_1_pH23_pH_3,zab2,zab3

      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)
     &                    +za(j1,j4)*zb(j4,j5)

      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s4(j1,j2,j3,j4)=s(j1,j2)+s(j1,j3)+s(j1,j4)
     &               +s(j2,j3)+s(j2,j4)+s(j3,j4)
C--- end statement functions

      qsq=
     & +s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j1,j5)
     &          +s(j2,j3)+s(j2,j4)+s(j2,j5)
     &                   +s(j3,j4)+s(j3,j5)
     &                            +s(j4,j5)
      sH1 = s4(j2,j3,j4,j5)
      sH2 = s4(j1,j3,j4,j5)
      sH3 = s4(j1,j2,j4,j5)
      s125 = s3(j1,j2,j5)
      s145 = s3(j1,j4,j5)
      s234 = s3(j2,j3,j4)
      s345 = s3(j3,j4,j5)

      zab_2_pH3_4 = -zab2(j2,j1,j5,j4)
      zab_5_pH4_3 = -zab2(j5,j1,j2,j3)
      zab_2_34_5 = zab2(j2,j3,j4,j5)
      zab_4_23_1 = zab2(j4,j2,j3,j1)
      zab_2_pH_1 = -zab3(j2,j3,j4,j5,j1)
      zab_2_pH_3 = -zab3(j2,j1,j4,j5,j3)
      zbb_1_pH_45_3 = -s345 * zb(j1,j3) - zb(j1,j2)*zab2(j2,j4,j5,j3)
      zbb_1_pH23_pH_3 = s145 * zb(j1,j3) + zb(j2,j3)*zab2(j2,j4,j5,j1)

      nAmmppp=
     & -za(j4,j5)**3/(za(j1,j2)*za(j1,j5)*za(j2,j3)*za(j3,j4))
     & -zab_5_pH4_3**3
     & /(s125*za(j1,j2)*za(j1,j5)*zab_2_pH3_4*zb(j4,j3))
     & -zab_4_23_1**3/(s234*za(j2,j3)*za(j3,j4)*zab_2_34_5*zb(j5,j1))
     & -qsq**2*zb(j3,j1)**4
     & /(sH2*zab_2_pH_1*zab_2_pH_3*zb(j4,j3)*zb(j5,j1)*zb(j5,j4))
     & -zbb_1_pH_45_3**3
     & /(sH1*s345*zab_2_pH_1*zab_2_34_5*zb(j4,j3)*zb(j5,j4))
     & +zbb_1_pH23_pH_3**3
     & /(sH3*s145*zab_2_pH_3*zab_2_pH3_4*zb(j5,j1)*zb(j5,j4))

      return
      end
