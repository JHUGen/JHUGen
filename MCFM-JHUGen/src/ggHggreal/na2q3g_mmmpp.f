      function na2q3g_mmmpp(j2,j3,j4,j5,j1,zb,za)
c--- Calcculation of the amplitudes using results of S. Badger
c--- as adapted from routines written by F. Caola      
C Arguments in call are strange because this
C has been created from A0Hqbqggg_mp_ppm
      
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::na2q3g_mmmpp
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::j1,j2,j3,j4,j5
      real(dp)::qsq,s3,s4,
     & sH1,sH3,sH4,s345,s234,s123,s125,s145
      complex(dp)::zab_2_pH1_5,zab_2_pH5_1,zab_2_15_4,zab_4_pH5_1,
     & zab_4_23_1,zab_5_pH_1,zab_5_pH_4,zab_5_12_3,
     & zbb_1_pH_45_3,zbb_1_pH_23_1,zbb_1_45_pH_3,
     & zbb_1_23_pH_4,zbb_3_pH_45_1,zbb_3_pH_15_4,zbb_3_pH12_pH_1,
     & zab2,zab3

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


      sH1 = s4(j2,j3,j4,j5)
      sH3 = s4(j1,j2,j4,j5)
      sH4 = s4(j1,j2,j3,j5)
      s345 = s3(j3,j4,j5)
      s234 = s3(j2,j3,j4)
      s123 = s3(j1,j2,j3)
      s125 = s3(j1,j2,j5)
      s145 = s3(j1,j4,j5)

      zab_2_pH1_5 = -zab2(j2,j3,j4,j5)
      zab_2_pH5_1 = -zab2(j2,j3,j4,j1)
      zab_4_pH5_1 = -zab2(j4,j2,j3,j1)
      zab_2_15_4 = zab2(j2,j1,j5,j4)
      zab_5_pH_1 = -zab3(j5,j2,j3,j4,j1)
      zab_5_pH_4 = -zab3(j5,j1,j2,j3,j4)
      zab_4_23_1 = zab2(j4,j2,j3,j1)
      zab_5_12_3 = zab2(j5,j1,j2,j3)

      zbb_1_pH_45_3=-zb(j1,j3)*s345-zb(j1,j2)*zab2(j2,j4,j5,j3)
      zbb_1_pH_23_1=
     & zb(j4,j1)*zab2(j4,j2,j3,j1)+zb(j5,j1)*zab2(j5,j2,j3,j1)
      zbb_1_45_pH_3=-zb(j1,j3)*s145-zb(j2,j3)*zab2(j2,j4,j5,j1)
      zbb_1_23_pH_4=-zb(j1,j4)*s123-zb(j5,j4)*zab2(j5,j2,j3,j1)
      zbb_3_pH_45_1=-zbb_1_45_pH_3
      zbb_3_pH_15_4=-zb(j3,j4)*s145-zb(j3,j2)*zab2(j2,j1,j5,j4)
      zbb_3_pH12_pH_1=zbb_1_pH_45_3

      na2q3g_mmmpp=
     & -za(j1,j5)**2*za(j2,j5)/(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5))
     & -za(j1,j5)**2*za(j2,j5)*zb(j4,j3)**3
     & /(s125*za(j1,j2)*zab_2_15_4*zab_5_12_3)
     & -s234**2*zab_2_pH5_1
     & /(za(j2,j3)*za(j3,j4)*zab_2_pH1_5*zab_4_pH5_1*zb(j5,j1))
     & +zab_5_pH_1**3*zb(j3,j1)*zb(j3,j2)**2
     & /(za(j4,j5)*zab_4_23_1*zb(j2,j1)*zbb_1_pH_45_3*zbb_1_pH_23_1)
     & -zab_5_pH_4**3*zb(j3,j1)*zb(j3,j2)**2
     & /(sH4*s123*zab_5_12_3*zb(j2,j1)*zbb_1_23_pH_4)
     & -qsq**2*zb(j3,j1)*zb(j3,j2)**2*zb(j4,j1)**4
     & /(zb(j2,j1)*zb(j5,j1)*zb(j5,j4)*zbb_1_pH_23_1*zbb_1_45_pH_3
     & *zbb_1_23_pH_4)
     & +zb(j4,j1)*zbb_3_pH_15_4**3
     & /(sH3*s145*zab_2_15_4*zb(j5,j1)*zb(j5,j4)*zbb_3_pH_45_1)
     & +sH1**2*zb(j4,j3)**3
     & /(s345*zab_2_pH1_5*zb(j5,j4)*zbb_3_pH12_pH_1)

      return
      end
