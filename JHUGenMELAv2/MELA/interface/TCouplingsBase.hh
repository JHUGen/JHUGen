#ifndef TCOUPLINGSBASE_HH
#define TCOUPLINGSBASE_HH


//---------------------------------
// Coupling array sizes
//---------------------------------
namespace{
  enum{
    gHIGGS_KAPPA,
    gHIGGS_KAPPA_TILDE,

    SIZE_HQQ
  };
  enum{
    gHIGGS_GG_2,
    gHIGGS_GG_3,
    gHIGGS_GG_4,

    SIZE_HGG
  };
  enum{
    gHIGGS_VV_1,
    gHIGGS_VV_2,
    gHIGGS_VV_3,
    gHIGGS_VV_4,

    gHIGGS_ZA_2,
    gHIGGS_ZA_3,
    gHIGGS_ZA_4,

    gHIGGS_AA_2,
    gHIGGS_AA_3,
    gHIGGS_AA_4,

    gHIGGS_VV_1_PRIME,
    gHIGGS_VV_1_PRIME2,
    gHIGGS_VV_1_PRIME3,
    gHIGGS_VV_1_PRIME4,
    gHIGGS_VV_1_PRIME5,

    gHIGGS_VV_2_PRIME,
    gHIGGS_VV_2_PRIME2,
    gHIGGS_VV_2_PRIME3,
    gHIGGS_VV_2_PRIME4,
    gHIGGS_VV_2_PRIME5,

    gHIGGS_VV_3_PRIME,
    gHIGGS_VV_3_PRIME2,
    gHIGGS_VV_3_PRIME3,
    gHIGGS_VV_3_PRIME4,
    gHIGGS_VV_3_PRIME5,

    gHIGGS_VV_4_PRIME,
    gHIGGS_VV_4_PRIME2,
    gHIGGS_VV_4_PRIME3,
    gHIGGS_VV_4_PRIME4,
    gHIGGS_VV_4_PRIME5,

    gHIGGS_ZA_1_PRIME2,

    gHIGGS_VV_1_PRIME6,
    gHIGGS_VV_1_PRIME7,
    gHIGGS_VV_2_PRIME6,
    gHIGGS_VV_2_PRIME7,
    gHIGGS_VV_3_PRIME6,
    gHIGGS_VV_3_PRIME7,
    gHIGGS_VV_4_PRIME6,
    gHIGGS_VV_4_PRIME7,

    SIZE_HVV
  };
  enum{
    LambdaHIGGS_QSQ_VV_1 = 0,
    LambdaHIGGS_QSQ_VV_2 = 1,
    LambdaHIGGS_QSQ_VV_3 = 2,
    LambdaHIGGS_QSQ_VV_4 = 3,

    SIZE_HVV_LAMBDAQSQ = 4
  };
  enum{
    cLambdaHIGGS_VV_QSQ1 = 0,
    cLambdaHIGGS_VV_QSQ2 = 1,
    cLambdaHIGGS_VV_QSQ12 = 2,

    SIZE_HVV_CQSQ = 3
  };
  enum{
    gZPRIME_QQ_LEFT,
    gZPRIME_QQ_RIGHT,

    SIZE_ZQQ
  };
  enum{
    gZPRIME_VV_1,
    gZPRIME_VV_2,

    SIZE_ZVV
  };
  enum{
    gGRAVITON_QQ_LEFT,
    gGRAVITON_QQ_RIGHT,

    SIZE_GQQ
  };
  enum{
    gGRAVITON_GG_1,
    gGRAVITON_GG_2,
    gGRAVITON_GG_3,
    gGRAVITON_GG_4,
    gGRAVITON_GG_5,

    SIZE_GGG
  };
  enum{
    gGRAVITON_VV_1,
    gGRAVITON_VV_2,
    gGRAVITON_VV_3,
    gGRAVITON_VV_4,
    gGRAVITON_VV_5,
    gGRAVITON_VV_6,
    gGRAVITON_VV_7,
    gGRAVITON_VV_8,
    gGRAVITON_VV_9,
    gGRAVITON_VV_10,

    SIZE_GVV
  };
}


#endif
