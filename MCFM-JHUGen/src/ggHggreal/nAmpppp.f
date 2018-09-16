      function nAmpppp(j5,j1,j2,j3,j4)
c--- Calculation of the amplitudes using results of S. Badger
c--- as adapted from routines written by F. Caola
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::nAmpppp
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::j1,j2,j3,j4,j5
      complex(dp)::zaa_5_pH_12_5, zaa_2_pH_34_5, zaa_2_pH_pH12_5,
     & zaa_3_pH_12_5,zab_5_pH_1, zab_2_pH_1, zab_5_pH_4, zab_3_pH_4,
     & zab_5_pH2_1, zab_5_pH3_4, zab_2_pH3_4, zab_3_pH2_1,
     & zab2, zab3
      real(dp)::sH1, sH2, sH3, sH4, s125, s145, s345
      real(dp)::s3, s4
      real(dp)::qsq
C--- begin statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=
     & za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)
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
      sH4 = s4(j1,j2,j3,j5)
      s125 = s3(j1,j2,j5)
      s145 = s3(j1,j4,j5)
      s345 = s3(j3,j4,j5)

      zaa_5_pH_12_5=-za(j1,j5)*zab2(j5,j3,j4,j1)
     & -za(j2,j5)*zab2(j5,j3,j4,j2)
      zaa_2_pH_34_5 = -za(j2,j5)*s345 - za(j2,j1)*zab2(j5,j3,j4,j1)
      zaa_2_pH_pH12_5 = - zaa_2_pH_34_5
      zaa_3_pH_12_5 = -za(j3,j5)*s125 - za(j3,j4)*zab2(j5,j1,j2,j4)

      zab_5_pH_1 = -zab3(j5,j2,j3,j4,j1)
      zab_2_pH_1 = -zab3(j2,j3,j4,j5,j1)
      zab_5_pH_4 = -zab3(j5,j1,j2,j3,j4)
      zab_3_pH_4 = -zab3(j3,j1,j2,j5,j4)
      zab_5_pH2_1 = -zab2(j5,j3,j4,j1)
      zab_5_pH3_4 = -zab2(j5,j1,j2,j4)
      zab_2_pH3_4 = -zab2(j2,j1,j5,j4)
      zab_3_pH2_1 = -zab2(j3,j4,j5,j1)

      nAmpppp = -zaa_5_pH_12_5**3/(za(j1,j2)*za(j3,j4)*za(j4,j5)
     & *za(j1,j5)*
     &   zaa_2_pH_34_5 * zaa_3_pH_12_5)
     &   + zab_5_pH_1**3/(sH1*za(j2,j3)*za(j3,j4)*za(j4,j5)*zab_2_pH_1)
     &   - zab_5_pH_4**3/(sH4*za(j1,j2)*za(j2,j3)*za(j1,j5)*zab_3_pH_4)
     &   - qsq**2 * zab_5_pH2_1**4/(sH2*s345*za(j3,j4)*za(j4,j5)*
     &   zaa_2_pH_pH12_5*zab_2_pH_1*zab_3_pH2_1)
     &   + qsq**2 * zab_5_pH3_4**4/(sH3*s125*za(j1,j2)*za(j1,j5)*
     &   zaa_3_pH_12_5*zab_2_pH3_4*zab_3_pH_4)
     &   -qsq**2 * zb(j4,j1)**4/(s145*za(j2,j3)*zab_2_pH3_4*
     &   zab_3_pH2_1*zb(j5,j1)*zb(j5,j4))

      return
      end

