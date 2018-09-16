      function na2q3g_mpppp(j2,j1,j3,j4,j5,za,zb)
c--- Calculation of the amplitudes using results of S. Badger
c--- as adapted from routines written by F. Caola      
C     A0Hqbqggg_mp_ppp
      
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::na2q3g_mpppp
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::j1,j2,j3,j4,j5
      complex(dp):: zaa_1_pH_34_2,zaa_2_15_pH_4,zaa_5_pH_34_2
      complex(dp):: zab_1_pH_3,zab_4_pH_3,zab_5_pH_1,zab_2_pH5_1,
     & zab_4_pH5_1,zab_1_25_3,zab_2_15_3,zab_5_12_3,
     & zab_4_23_1,zab2, zab3
      real(dp)::qsq,sH1,sH3,sH4,sH5,s123,s125,s234,s3,s4
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
      sH5 = s4(j1,j2,j3,j4)
      s123 = s3(j1,j2,j3)
      s125 = s3(j1,j2,j5)
      s234 = s3(j2,j3,j4)

      zaa_1_pH_34_2 = -za(j1,j2)*s234 - za(j1,j5)*zab2(j2,j3,j4,j5)
      zaa_2_15_pH_4 = -za(j2,j4)*s125 - za(j3,j4)*zab2(j2,j1,j5,j3)
      zaa_5_pH_34_2 = -za(j5,j2)*s234 - za(j5,j1)*zab2(j2,j3,j4,j1)

      zab_1_pH_3 = -zab3(j1,j2,j4,j5,j3)
      zab_4_pH_3 = -zab3(j4,j1,j2,j5,j3)
      zab_5_pH_1 = -zab3(j5,j2,j3,j4,j1)
      zab_2_pH5_1 = -zab2(j2,j3,j4,j1)
      zab_4_pH5_1 = -zab2(j4,j2,j3,j1)

      zab_1_25_3 = zab2(j1,j2,j5,j3)
      zab_2_15_3 = zab2(j2,j1,j5,j3)
      zab_5_12_3 = zab2(j5,j1,j2,j3)
      zab_4_23_1 = -zab_4_pH5_1

      na2q3g_mpppp=zaa_1_pH_34_2**3
     & /(za(j1,j2)*za(j1,j5)*za(j2,j3)*za(j3,j4)
     & *zaa_2_15_pH_4*zaa_5_pH_34_2)
     & -zab_1_pH_3**3
     & /(sH3*za(j1,j2)*za(j1,j5)*za(j4,j5)*zab_4_pH_3)
     & +sH1**2/(za(j2,j3)*za(j3,j4)*za(j4,j5)*zab_5_pH_1)
     & -qsq**2*s234**2*zab_2_pH5_1
     & /(sH5*za(j2,j3)*za(j3,j4)*zaa_5_pH_34_2*zab_4_pH5_1*zab_5_pH_1)
     & -qsq**2*zab_1_25_3**3*zab_2_15_3
     & /(sH4*s125*za(j1,j2)*za(j1,j5)*zaa_2_15_pH_4
     & *zab_4_pH_3*zab_5_12_3)
     & -qsq**2*zb(j3,j1)*zb(j3,j2)**2
     & /(s123*za(j4,j5)*zab_4_23_1*zab_5_12_3*zb(j2,j1))
      return
      end

