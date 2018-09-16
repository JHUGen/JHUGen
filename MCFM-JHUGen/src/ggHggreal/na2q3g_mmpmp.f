      function na2q3g_mmpmp(j2,j3,j4,j5,j1,zb,za)
c--- Calculation of the amplitudes using results of S. Badger
c--- as adapted from routines written by F. Caola      
C Arguments in call are strange because this 
C has been created from A0Hqbqggg_mp_pmp
      
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::na2q3g_mmpmp
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::j1,j2,j3,j4,j5
      real(dp)::qsq,s123,s125,s145,s234,s345,sH1,sH3,sH5,s3,s4
      complex(dp)::zab_2_pH1_5,zab_4_pH1_5,zab_4_pH5_1,zab_4_pH_5,
     & zab_4_23_1,zab_4_pH_3,zab_1_25_3,zab_2_15_3,zab_2_15_4,
     & zab_5_12_3,zbb_3_pH_14_5,zbb_3_pH_45_1,zbb_3_12_pH_5,
     & zbb_1_45_pH_3,zbb_3_45_pH_3,zbb_1_pH_45_3,zbb_3_pH12_pH_1,
     & zab2, zab3

c--- begin statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=
     & za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)

      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s4(j1,j2,j3,j4)=s(j1,j2)+s(j1,j3)+s(j1,j4)
     &               +s(j2,j3)+s(j2,j4)+s(j3,j4)
c--- end statement functions

      qsq=
     & +s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j1,j5)
     &          +s(j2,j3)+s(j2,j4)+s(j2,j5)
     &                   +s(j3,j4)+s(j3,j5)
     &                            +s(j4,j5)

      s123 = s3(j1,j2,j3)
      s125 = s3(j1,j2,j5)
      s145 = s3(j1,j4,j5)
      s234 = s3(j2,j3,j4)
      s345 = s3(j3,j4,j5)
      sH1 = s4(j2,j3,j4,j5)
      sH3 = s4(j1,j2,j4,j5)
      sH5 = s4(j1,j2,j3,j4)

      zab_2_pH1_5 = -zab2(j2,j3,j4,j5)
      zab_4_pH1_5 = -zab2(j4,j2,j3,j5)
      zab_4_pH5_1 = -zab2(j4,j2,j3,j1)
      zab_4_pH_5 = -zab3(j4,j1,j2,j3,j5)
      zab_4_pH_3 = -zab3(j4,j1,j2,j5,j3)

      zab_1_25_3 = zab2(j1,j2,j5,j3)
      zab_4_23_1 = zab2(j4,j2,j3,j1)
      zab_2_15_3 = zab2(j2,j1,j5,j3)
      zab_2_15_4 = zab2(j2,j1,j5,j4)
      zab_5_12_3 = zab2(j5,j1,j2,j3)

      zbb_3_pH_14_5 = -zb(j3,j5)*s145 - zb(j3,j2)*zab2(j2,j1,j4,j5)
      zbb_3_pH_45_1 = -zb(j3,j1)*s145 - zb(j3,j2)*zab2(j2,j4,j5,j1)
      zbb_3_12_pH_5 = -zb(j3,j5)*s123 - zb(j4,j5)*zab2(j4,j1,j2,j3)
      zbb_1_45_pH_3 = -zb(j1,j3)*s145 - zb(j2,j3)*zab2(j2,j4,j5,j1)
      zbb_3_45_pH_3 =
     &  zb(j4,j3)*zab2(j4,j1,j2,j3) + zb(j5,j3)*zab2(j5,j1,j2,j3)
      zbb_1_pH_45_3 = -zb(j1,j3)*s345 - zb(j1,j2)*zab2(j2,j4,j5,j3)
      zbb_3_pH12_pH_1 = zb(j3,j1)*s345 + zb(j2,j1)*zab2(j2,j4,j5,j3)

      na2q3g_mmpmp= -za(j1,j4)**3*za(j2,j4)
     & /(za(j1,j2)*za(j1,j5)*za(j2,j3)*za(j3,j4)*za(j4,j5))
     & -za(j2,j4)*zab_4_pH1_5**3
     & /(s234*za(j2,j3)*za(j3,j4)*zab_2_pH1_5*zab_4_pH5_1)
     & -zab_1_25_3**3*zab_2_15_3
     & /(s125*za(j1,j2)*za(j1,j5)*zab_2_15_4*zab_5_12_3*zb(j4,j3))
     & +zbb_3_pH_14_5**3
     & /(sH3*s145*zab_2_15_4*zb(j5,j4)*zbb_3_pH_45_1)
     & +zab_4_pH_5**3*zb(j3,j1)*zb(j3,j2)**2 
     & /(sH5*s123*zab_4_23_1*zb(j2,j1)*zbb_3_12_pH_5)
     & -zab_4_pH_3**3*zb(j3,j1)*zb(j3,j2)**2
     & /(za(j4,j5)*zab_5_12_3*zb(j2,j1)*zbb_1_45_pH_3*zbb_3_45_pH_3)
     & +qsq**2*zb(j3,j1)*zb(j3,j2)**2*zb(j5,j3)**4
     & /(zb(j2,j1)*zb(j4,j3)*zb(j5,j4)*zbb_1_pH_45_3
     & *zbb_3_12_pH_5*zbb_3_45_pH_3)
     & +sH1**2*zb(j5,j3)**4
     & /(s345*zab_2_pH1_5*zb(j4,j3)*zb(j5,j4)*zbb_3_pH12_pH_1)

      return
      end

