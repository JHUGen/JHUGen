#include "RooSpinZero.h"


using namespace std;
using namespace MELAStreamHelpers;


RooSpinZero::RooSpinZero() : RooSpin(){}
RooSpinZero::RooSpinZero(
  const char* name, const char* title,
  modelMeasurables const& _measurables,
  modelParameters const& _parameters,
  modelCouplings const& _couplings,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2,
  TVar::VerbosityLevel verbosity_
) : RooSpin(
  name, title,
  _measurables, _parameters,
  _Vdecay1, _Vdecay2,
  verbosity_
  ),

  g1Val("g1Val", "g1Val", this, (RooAbsReal&)*(_couplings.g1List[0][0])),
  g2Val("g2Val", "g2Val", this, (RooAbsReal&)*(_couplings.g2List[0][0])),
  g3Val("g3Val", "g3Val", this, (RooAbsReal&)*(_couplings.g3List[0][0])),
  g4Val("g4Val", "g4Val", this, (RooAbsReal&)*(_couplings.g4List[0][0])),

  g1_primeVal("g1_primeVal", "g1_primeVal", this, (RooAbsReal&)*(_couplings.g1List[1][0])),
  g2_primeVal("g2_primeVal", "g2_primeVal", this, (RooAbsReal&)*(_couplings.g2List[1][0])),
  g3_primeVal("g3_primeVal", "g3_primeVal", this, (RooAbsReal&)*(_couplings.g3List[1][0])),
  g4_primeVal("g4_primeVal", "g4_primeVal", this, (RooAbsReal&)*(_couplings.g4List[1][0])),

  g1_prime2Val("g1_prime2Val", "g1_prime2Val", this, (RooAbsReal&)*(_couplings.g1List[2][0])),
  g2_prime2Val("g2_prime2Val", "g2_prime2Val", this, (RooAbsReal&)*(_couplings.g2List[2][0])),
  g3_prime2Val("g3_prime2Val", "g3_prime2Val", this, (RooAbsReal&)*(_couplings.g3List[2][0])),
  g4_prime2Val("g4_prime2Val", "g4_prime2Val", this, (RooAbsReal&)*(_couplings.g4List[2][0])),

  g1_prime3Val("g1_prime3Val", "g1_prime3Val", this, (RooAbsReal&)*(_couplings.g1List[3][0])),
  g2_prime3Val("g2_prime3Val", "g2_prime3Val", this, (RooAbsReal&)*(_couplings.g2List[3][0])),
  g3_prime3Val("g3_prime3Val", "g3_prime3Val", this, (RooAbsReal&)*(_couplings.g3List[3][0])),
  g4_prime3Val("g4_prime3Val", "g4_prime3Val", this, (RooAbsReal&)*(_couplings.g4List[3][0])),

  g1_prime4Val("g1_prime4Val", "g1_prime4Val", this, (RooAbsReal&)*(_couplings.g1List[4][0])),
  g2_prime4Val("g2_prime4Val", "g2_prime4Val", this, (RooAbsReal&)*(_couplings.g2List[4][0])),
  g3_prime4Val("g3_prime4Val", "g3_prime4Val", this, (RooAbsReal&)*(_couplings.g3List[4][0])),
  g4_prime4Val("g4_prime4Val", "g4_prime4Val", this, (RooAbsReal&)*(_couplings.g4List[4][0])),

  g1_prime5Val("g1_prime5Val", "g1_prime5Val", this, (RooAbsReal&)*(_couplings.g1List[5][0])),
  g2_prime5Val("g2_prime5Val", "g2_prime5Val", this, (RooAbsReal&)*(_couplings.g2List[5][0])),
  g3_prime5Val("g3_prime5Val", "g3_prime5Val", this, (RooAbsReal&)*(_couplings.g3List[5][0])),
  g4_prime5Val("g4_prime5Val", "g4_prime5Val", this, (RooAbsReal&)*(_couplings.g4List[5][0])),

  g1_prime6Val("g1_prime6Val", "g1_prime6Val", this, (RooAbsReal&)*(_couplings.g1List[6][0])),
  g2_prime6Val("g2_prime6Val", "g2_prime6Val", this, (RooAbsReal&)*(_couplings.g2List[6][0])),
  g3_prime6Val("g3_prime6Val", "g3_prime6Val", this, (RooAbsReal&)*(_couplings.g3List[6][0])),
  g4_prime6Val("g4_prime6Val", "g4_prime6Val", this, (RooAbsReal&)*(_couplings.g4List[6][0])),

  g1_prime7Val("g1_prime7Val", "g1_prime7Val", this, (RooAbsReal&)*(_couplings.g1List[7][0])),
  g2_prime7Val("g2_prime7Val", "g2_prime7Val", this, (RooAbsReal&)*(_couplings.g2List[7][0])),
  g3_prime7Val("g3_prime7Val", "g3_prime7Val", this, (RooAbsReal&)*(_couplings.g3List[7][0])),
  g4_prime7Val("g4_prime7Val", "g4_prime7Val", this, (RooAbsReal&)*(_couplings.g4List[7][0])),

  gzgs1_prime2Val("gzgs1_prime2Val", "gzgs1_prime2Val", this, (RooAbsReal&)*(_couplings.gzgs1List[0][0])), // Special case!
  gzgs2Val("gzgs2Val", "gzgs2Val", this, (RooAbsReal&)*(_couplings.gzgs2List[0][0])),
  gzgs3Val("gzgs3Val", "gzgs3Val", this, (RooAbsReal&)*(_couplings.gzgs3List[0][0])),
  gzgs4Val("gzgs4Val", "gzgs4Val", this, (RooAbsReal&)*(_couplings.gzgs4List[0][0])),
  ggsgs2Val("ggsgs2Val", "ggsgs2Val", this, (RooAbsReal&)*(_couplings.ggsgs2List[0][0])),
  ggsgs3Val("ggsgs3Val", "ggsgs3Val", this, (RooAbsReal&)*(_couplings.ggsgs3List[0][0])),
  ggsgs4Val("ggsgs4Val", "ggsgs4Val", this, (RooAbsReal&)*(_couplings.ggsgs4List[0][0])),

  g1ValIm("g1ValIm", "g1ValIm", this, (RooAbsReal&)*(_couplings.g1List[0][1])),
  g2ValIm("g2ValIm", "g2ValIm", this, (RooAbsReal&)*(_couplings.g2List[0][1])),
  g3ValIm("g3ValIm", "g3ValIm", this, (RooAbsReal&)*(_couplings.g3List[0][1])),
  g4ValIm("g4ValIm", "g4ValIm", this, (RooAbsReal&)*(_couplings.g4List[0][1])),

  g1_primeValIm("g1_primeValIm", "g1_primeValIm", this, (RooAbsReal&)*(_couplings.g1List[1][1])),
  g2_primeValIm("g2_primeValIm", "g2_primeValIm", this, (RooAbsReal&)*(_couplings.g2List[1][1])),
  g3_primeValIm("g3_primeValIm", "g3_primeValIm", this, (RooAbsReal&)*(_couplings.g3List[1][1])),
  g4_primeValIm("g4_primeValIm", "g4_primeValIm", this, (RooAbsReal&)*(_couplings.g4List[1][1])),

  g1_prime2ValIm("g1_prime2ValIm", "g1_prime2ValIm", this, (RooAbsReal&)*(_couplings.g1List[2][1])),
  g2_prime2ValIm("g2_prime2ValIm", "g2_prime2ValIm", this, (RooAbsReal&)*(_couplings.g2List[2][1])),
  g3_prime2ValIm("g3_prime2ValIm", "g3_prime2ValIm", this, (RooAbsReal&)*(_couplings.g3List[2][1])),
  g4_prime2ValIm("g4_prime2ValIm", "g4_prime2ValIm", this, (RooAbsReal&)*(_couplings.g4List[2][1])),

  g1_prime3ValIm("g1_prime3ValIm", "g1_prime3ValIm", this, (RooAbsReal&)*(_couplings.g1List[3][1])),
  g2_prime3ValIm("g2_prime3ValIm", "g2_prime3ValIm", this, (RooAbsReal&)*(_couplings.g2List[3][1])),
  g3_prime3ValIm("g3_prime3ValIm", "g3_prime3ValIm", this, (RooAbsReal&)*(_couplings.g3List[3][1])),
  g4_prime3ValIm("g4_prime3ValIm", "g4_prime3ValIm", this, (RooAbsReal&)*(_couplings.g4List[3][1])),

  g1_prime4ValIm("g1_prime4ValIm", "g1_prime4ValIm", this, (RooAbsReal&)*(_couplings.g1List[4][1])),
  g2_prime4ValIm("g2_prime4ValIm", "g2_prime4ValIm", this, (RooAbsReal&)*(_couplings.g2List[4][1])),
  g3_prime4ValIm("g3_prime4ValIm", "g3_prime4ValIm", this, (RooAbsReal&)*(_couplings.g3List[4][1])),
  g4_prime4ValIm("g4_prime4ValIm", "g4_prime4ValIm", this, (RooAbsReal&)*(_couplings.g4List[4][1])),

  g1_prime5ValIm("g1_prime5ValIm", "g1_prime5ValIm", this, (RooAbsReal&)*(_couplings.g1List[5][1])),
  g2_prime5ValIm("g2_prime5ValIm", "g2_prime5ValIm", this, (RooAbsReal&)*(_couplings.g2List[5][1])),
  g3_prime5ValIm("g3_prime5ValIm", "g3_prime5ValIm", this, (RooAbsReal&)*(_couplings.g3List[5][1])),
  g4_prime5ValIm("g4_prime5ValIm", "g4_prime5ValIm", this, (RooAbsReal&)*(_couplings.g4List[5][1])),

  g1_prime6ValIm("g1_prime6ValIm", "g1_prime6ValIm", this, (RooAbsReal&)*(_couplings.g1List[6][1])),
  g2_prime6ValIm("g2_prime6ValIm", "g2_prime6ValIm", this, (RooAbsReal&)*(_couplings.g2List[6][1])),
  g3_prime6ValIm("g3_prime6ValIm", "g3_prime6ValIm", this, (RooAbsReal&)*(_couplings.g3List[6][1])),
  g4_prime6ValIm("g4_prime6ValIm", "g4_prime6ValIm", this, (RooAbsReal&)*(_couplings.g4List[6][1])),

  g1_prime7ValIm("g1_prime7ValIm", "g1_prime7ValIm", this, (RooAbsReal&)*(_couplings.g1List[7][1])),
  g2_prime7ValIm("g2_prime7ValIm", "g2_prime7ValIm", this, (RooAbsReal&)*(_couplings.g2List[7][1])),
  g3_prime7ValIm("g3_prime7ValIm", "g3_prime7ValIm", this, (RooAbsReal&)*(_couplings.g3List[7][1])),
  g4_prime7ValIm("g4_prime7ValIm", "g4_prime7ValIm", this, (RooAbsReal&)*(_couplings.g4List[7][1])),

  gzgs1_prime2ValIm("gzgs1_prime2ValIm", "gzgs1_prime2ValIm", this, (RooAbsReal&)*(_couplings.gzgs1List[0][1])), // Special case!
  gzgs2ValIm("gzgs2ValIm", "gzgs2ValIm", this, (RooAbsReal&)*(_couplings.gzgs2List[0][1])),
  gzgs3ValIm("gzgs3ValIm", "gzgs3ValIm", this, (RooAbsReal&)*(_couplings.gzgs3List[0][1])),
  gzgs4ValIm("gzgs4ValIm", "gzgs4ValIm", this, (RooAbsReal&)*(_couplings.gzgs4List[0][1])),
  ggsgs2ValIm("ggsgs2ValIm", "ggsgs2ValIm", this, (RooAbsReal&)*(_couplings.ggsgs2List[0][1])),
  ggsgs3ValIm("ggsgs3ValIm", "ggsgs3ValIm", this, (RooAbsReal&)*(_couplings.ggsgs3List[0][1])),
  ggsgs4ValIm("ggsgs4ValIm", "ggsgs4ValIm", this, (RooAbsReal&)*(_couplings.ggsgs4List[0][1])),

  Lambda("Lambda", "Lambda", this, (RooAbsReal&)*(_couplings.Lambda)),
  Lambda_zgs1("Lambda_zgs1", "Lambda_zgs1", this, (RooAbsReal&)*(_couplings.Lambda_zgs1)),
  Lambda_z1("Lambda_z1", "Lambda_z1", this, (RooAbsReal&)*(_couplings.Lambda_z1)),
  Lambda_z2("Lambda_z2", "Lambda_z2", this, (RooAbsReal&)*(_couplings.Lambda_z2)),
  Lambda_z3("Lambda_z3", "Lambda_z3", this, (RooAbsReal&)*(_couplings.Lambda_z3)),
  Lambda_z4("Lambda_z4", "Lambda_z4", this, (RooAbsReal&)*(_couplings.Lambda_z4)),
  Lambda_Q("Lambda_Q", "Lambda_Q", this, (RooAbsReal&)*(_couplings.Lambda_Q)),

  Lambda_z11("Lambda_z11", "Lambda_z11", this, (RooAbsReal&)*(_couplings.Lambda_z1qsq[0])),
  Lambda_z21("Lambda_z21", "Lambda_z21", this, (RooAbsReal&)*(_couplings.Lambda_z2qsq[0])),
  Lambda_z31("Lambda_z31", "Lambda_z31", this, (RooAbsReal&)*(_couplings.Lambda_z3qsq[0])),
  Lambda_z41("Lambda_z41", "Lambda_z41", this, (RooAbsReal&)*(_couplings.Lambda_z4qsq[0])),

  Lambda_z12("Lambda_z12", "Lambda_z12", this, (RooAbsReal&)*(_couplings.Lambda_z1qsq[1])),
  Lambda_z22("Lambda_z22", "Lambda_z22", this, (RooAbsReal&)*(_couplings.Lambda_z2qsq[1])),
  Lambda_z32("Lambda_z32", "Lambda_z32", this, (RooAbsReal&)*(_couplings.Lambda_z3qsq[1])),
  Lambda_z42("Lambda_z42", "Lambda_z42", this, (RooAbsReal&)*(_couplings.Lambda_z4qsq[1])),

  Lambda_z10("Lambda_z10", "Lambda_z10", this, (RooAbsReal&)*(_couplings.Lambda_z1qsq[2])),
  Lambda_z20("Lambda_z20", "Lambda_z20", this, (RooAbsReal&)*(_couplings.Lambda_z2qsq[2])),
  Lambda_z30("Lambda_z30", "Lambda_z30", this, (RooAbsReal&)*(_couplings.Lambda_z3qsq[2])),
  Lambda_z40("Lambda_z40", "Lambda_z40", this, (RooAbsReal&)*(_couplings.Lambda_z4qsq[2])),

  cz_q1sq("cz_q1sq", "cz_q1sq", this, (RooAbsReal&)*(_couplings.cLambda_qsq[0])),
  cz_q2sq("cz_q2sq", "cz_q2sq", this, (RooAbsReal&)*(_couplings.cLambda_qsq[1])),
  cz_q12sq("cz_q12sq", "cz_q12sq", this, (RooAbsReal&)*(_couplings.cLambda_qsq[2])),

  gvvp1Val("gvvp1Val", "gvvp1Val", this, (RooAbsReal&)*(_couplings.gvvp1List[0][0])),
  gvpvp1Val("gvpvp1Val", "gvpvp1Val", this, (RooAbsReal&)*(_couplings.gvpvp1List[0][0])),

  gvvp1ValIm("gvvp1ValIm", "gvvp1ValIm", this, (RooAbsReal&)*(_couplings.gvvp1List[0][1])),
  gvpvp1ValIm("gvpvp1ValIm", "gvpvp1ValIm", this, (RooAbsReal&)*(_couplings.gvpvp1List[0][1]))
{}


RooSpinZero::RooSpinZero(const RooSpinZero& other, const char* name) :
RooSpin(other, name),

g1Val("g1Val", this, other.g1Val),
g2Val("g2Val", this, other.g2Val),
g3Val("g3Val", this, other.g3Val),
g4Val("g4Val", this, other.g4Val),

g1_primeVal("g1_primeVal", this, other.g1_primeVal),
g2_primeVal("a2_primeVal", this, other.g2_primeVal),
g3_primeVal("g3_primeVal", this, other.g3_primeVal),
g4_primeVal("g4_primeVal", this, other.g4_primeVal),

g1_prime2Val("g1_prime2Val", this, other.g1_prime2Val),
g2_prime2Val("a2_prime2Val", this, other.g2_prime2Val),
g3_prime2Val("g3_prime2Val", this, other.g3_prime2Val),
g4_prime2Val("g4_prime2Val", this, other.g4_prime2Val),

g1_prime3Val("g1_prime3Val", this, other.g1_prime3Val),
g2_prime3Val("a2_prime3Val", this, other.g2_prime3Val),
g3_prime3Val("g3_prime3Val", this, other.g3_prime3Val),
g4_prime3Val("g4_prime3Val", this, other.g4_prime3Val),

g1_prime4Val("g1_prime4Val", this, other.g1_prime4Val),
g2_prime4Val("a2_prime4Val", this, other.g2_prime4Val),
g3_prime4Val("g3_prime4Val", this, other.g3_prime4Val),
g4_prime4Val("g4_prime4Val", this, other.g4_prime4Val),

g1_prime5Val("g1_prime5Val", this, other.g1_prime5Val),
g2_prime5Val("a2_prime5Val", this, other.g2_prime5Val),
g3_prime5Val("g3_prime5Val", this, other.g3_prime5Val),
g4_prime5Val("g4_prime5Val", this, other.g4_prime5Val),

g1_prime6Val("g1_prime6Val", this, other.g1_prime6Val),
g2_prime6Val("a2_prime6Val", this, other.g2_prime6Val),
g3_prime6Val("g3_prime6Val", this, other.g3_prime6Val),
g4_prime6Val("g4_prime6Val", this, other.g4_prime6Val),

g1_prime7Val("g1_prime7Val", this, other.g1_prime7Val),
g2_prime7Val("a2_prime7Val", this, other.g2_prime7Val),
g3_prime7Val("g3_prime7Val", this, other.g3_prime7Val),
g4_prime7Val("g4_prime7Val", this, other.g4_prime7Val),

gzgs1_prime2Val("gzgs1_prime2Val", this, other.gzgs1_prime2Val),
gzgs2Val("gzgs2Val", this, other.gzgs2Val),
gzgs3Val("gzgs3Val", this, other.gzgs3Val),
gzgs4Val("gzgs4Val", this, other.gzgs4Val),
ggsgs2Val("ggsgs2Val", this, other.ggsgs2Val),
ggsgs3Val("ggsgs3Val", this, other.ggsgs3Val),
ggsgs4Val("ggsgs4Val", this, other.ggsgs4Val),


g1ValIm("g1ValIm", this, other.g1ValIm),
g2ValIm("g2ValIm", this, other.g2ValIm),
g3ValIm("g3ValIm", this, other.g3ValIm),
g4ValIm("g4ValIm", this, other.g4ValIm),

g1_primeValIm("g1_primeValIm", this, other.g1_primeValIm),
g2_primeValIm("a2_primeValIm", this, other.g2_primeValIm),
g3_primeValIm("g3_primeValIm", this, other.g3_primeValIm),
g4_primeValIm("g4_primeValIm", this, other.g4_primeValIm),

g1_prime2ValIm("g1_prime2ValIm", this, other.g1_prime2ValIm),
g2_prime2ValIm("a2_prime2ValIm", this, other.g2_prime2ValIm),
g3_prime2ValIm("g3_prime2ValIm", this, other.g3_prime2ValIm),
g4_prime2ValIm("g4_prime2ValIm", this, other.g4_prime2ValIm),

g1_prime3ValIm("g1_prime3ValIm", this, other.g1_prime3ValIm),
g2_prime3ValIm("a2_prime3ValIm", this, other.g2_prime3ValIm),
g3_prime3ValIm("g3_prime3ValIm", this, other.g3_prime3ValIm),
g4_prime3ValIm("g4_prime3ValIm", this, other.g4_prime3ValIm),

g1_prime4ValIm("g1_prime4ValIm", this, other.g1_prime4ValIm),
g2_prime4ValIm("a2_prime4ValIm", this, other.g2_prime4ValIm),
g3_prime4ValIm("g3_prime4ValIm", this, other.g3_prime4ValIm),
g4_prime4ValIm("g4_prime4ValIm", this, other.g4_prime4ValIm),

g1_prime5ValIm("g1_prime5ValIm", this, other.g1_prime5ValIm),
g2_prime5ValIm("a2_prime5ValIm", this, other.g2_prime5ValIm),
g3_prime5ValIm("g3_prime5ValIm", this, other.g3_prime5ValIm),
g4_prime5ValIm("g4_prime5ValIm", this, other.g4_prime5ValIm),

g1_prime6ValIm("g1_prime6ValIm", this, other.g1_prime6ValIm),
g2_prime6ValIm("a2_prime6ValIm", this, other.g2_prime6ValIm),
g3_prime6ValIm("g3_prime6ValIm", this, other.g3_prime6ValIm),
g4_prime6ValIm("g4_prime6ValIm", this, other.g4_prime6ValIm),

g1_prime7ValIm("g1_prime7ValIm", this, other.g1_prime7ValIm),
g2_prime7ValIm("a2_prime7ValIm", this, other.g2_prime7ValIm),
g3_prime7ValIm("g3_prime7ValIm", this, other.g3_prime7ValIm),
g4_prime7ValIm("g4_prime7ValIm", this, other.g4_prime7ValIm),

gzgs1_prime2ValIm("gzgs1_prime2ValIm", this, other.gzgs1_prime2ValIm),
gzgs2ValIm("gzgs2ValIm", this, other.gzgs2ValIm),
gzgs3ValIm("gzgs3ValIm", this, other.gzgs3ValIm),
gzgs4ValIm("gzgs4ValIm", this, other.gzgs4ValIm),
ggsgs2ValIm("ggsgs2ValIm", this, other.ggsgs2ValIm),
ggsgs3ValIm("ggsgs3ValIm", this, other.ggsgs3ValIm),
ggsgs4ValIm("ggsgs4ValIm", this, other.ggsgs4ValIm),

Lambda("Lambda", this, other.Lambda),
Lambda_zgs1("Lambda_zgs1", this, other.Lambda_zgs1),
Lambda_z1("Lambda_z1", this, other.Lambda_z1),
Lambda_z2("Lambda_z2", this, other.Lambda_z2),
Lambda_z3("Lambda_z3", this, other.Lambda_z3),
Lambda_z4("Lambda_z4", this, other.Lambda_z4),
Lambda_Q("Lambda_Q", this, other.Lambda_Q),

Lambda_z11("Lambda_z11", this, other.Lambda_z11),
Lambda_z21("Lambda_z21", this, other.Lambda_z21),
Lambda_z31("Lambda_z31", this, other.Lambda_z31),
Lambda_z41("Lambda_z41", this, other.Lambda_z41),

Lambda_z12("Lambda_z12", this, other.Lambda_z12),
Lambda_z22("Lambda_z22", this, other.Lambda_z22),
Lambda_z32("Lambda_z32", this, other.Lambda_z32),
Lambda_z42("Lambda_z42", this, other.Lambda_z42),

Lambda_z10("Lambda_z10", this, other.Lambda_z10),
Lambda_z20("Lambda_z20", this, other.Lambda_z20),
Lambda_z30("Lambda_z30", this, other.Lambda_z30),
Lambda_z40("Lambda_z40", this, other.Lambda_z40),

cz_q1sq("cz_q1sq", this, other.cz_q1sq),
cz_q2sq("cz_q2sq", this, other.cz_q2sq),
cz_q12sq("cz_q12sq", this, other.cz_q12sq),

gvvp1Val("gvvp1Val", this, other.gvvp1Val),
gvpvp1Val("gvpvp1Val", this, other.gvpvp1Val),

gvvp1ValIm("gvvp1ValIm", this, other.gvvp1ValIm),
gvpvp1ValIm("gvpvp1ValIm", this, other.gvpvp1ValIm)
{}

void RooSpinZero::calculateAi(
  Double_t& a1Re, Double_t& a1Im, Double_t& a2Re, Double_t& a2Im, Double_t& a3Re, Double_t& a3Im,
  int VGammaVpmode1, int VGammaVpmode2
)const{
  Double_t mV;
  getMVGamV(&mV);

  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;

  Double_t s = (pow(m12, 2) - pow(m1_, 2) - pow(m2_, 2))/2.;
  if (m1_>m2_+m12 || m2_>m1_+m12) s = -s;

  Double_t g1_dyn=0;
  Double_t g2_dyn=0;
  Double_t g3_dyn=0;
  Double_t g4_dyn=0;
  Double_t g1_dynIm=0;
  Double_t g2_dynIm=0;
  Double_t g3_dynIm=0;
  Double_t g4_dynIm=0;

  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell){ // VV/VVp/VpVp
    if (VGammaVpmode1==0 && VGammaVpmode2==0){ // VV
      if (g1_primeVal!=0) g1_dyn += g1_primeVal * pow(Lambda_z1, 4)/(pow(Lambda_z1, 2) + pow(m1_, 2))/(pow(Lambda_z1, 2) + pow(m2_, 2));
      if (g1_prime2Val!=0) g1_dyn += g1_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z1, 2);
      if (g1_prime3Val!=0) g1_dyn += g1_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z1, 2);
      if (g1_prime4Val!=0) g1_dyn += g1_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
      if (g1_prime5Val!=0) g1_dyn += g1_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z1, 4);
      if (g1_prime6Val!=0) g1_dyn += g1_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z1, 4);
      if (g1_prime7Val!=0) g1_dyn += g1_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z1, 4);

      if (cz_q1sq!=0.) g1_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z11, 2));
      if (cz_q2sq!=0.) g1_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z12, 2));
      if (cz_q12sq!=0.) g1_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z10, 2));
      g1_dyn += g1Val;

      g2_dyn = g2Val;
      if (g2_primeVal!=0) g2_dyn += g2_primeVal * pow(Lambda_z2, 4)/(pow(Lambda_z2, 2) + pow(m1_, 2))/(pow(Lambda_z2, 2) + pow(m2_, 2));
      if (g2_prime2Val!=0) g2_dyn += g2_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z2, 2);
      if (g2_prime3Val!=0) g2_dyn += g2_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z2, 2);
      if (g2_prime4Val!=0) g2_dyn += g2_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
      if (g2_prime5Val!=0) g2_dyn += g2_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z2, 4);
      if (g2_prime6Val!=0) g2_dyn += g2_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z2, 4);
      if (g2_prime7Val!=0) g2_dyn += g2_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z2, 4);

      if (cz_q1sq!=0.) g2_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z21, 2));
      if (cz_q2sq!=0.) g2_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z22, 2));
      if (cz_q12sq!=0.) g2_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z20, 2));

      g3_dyn = g3Val;
      if (g3_primeVal!=0) g3_dyn += g3_primeVal * pow(Lambda_z3, 4)/(pow(Lambda_z3, 2) + pow(m1_, 2))/(pow(Lambda_z3, 2) + pow(m2_, 2));
      if (g3_prime2Val!=0) g3_dyn += g3_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z3, 2);
      if (g3_prime3Val!=0) g3_dyn += g3_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z3, 2);
      if (g3_prime4Val!=0) g3_dyn += g3_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
      if (g3_prime5Val!=0) g3_dyn += g3_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z3, 4);
      if (g3_prime6Val!=0) g3_dyn += g3_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z3, 4);
      if (g3_prime7Val!=0) g3_dyn += g3_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z3, 4);

      if (cz_q1sq!=0.) g3_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z31, 2));
      if (cz_q2sq!=0.) g3_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z32, 2));
      if (cz_q12sq!=0.) g3_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z30, 2));

      g4_dyn = g4Val;
      if (g4_primeVal!=0) g4_dyn += g4_primeVal * pow(Lambda_z4, 4)/(pow(Lambda_z4, 2) + pow(m1_, 2))/(pow(Lambda_z4, 2) + pow(m2_, 2));
      if (g4_prime2Val!=0) g4_dyn += g4_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z4, 2);
      if (g4_prime3Val!=0) g4_dyn += g4_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z4, 2);
      if (g4_prime4Val!=0) g4_dyn += g4_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
      if (g4_prime5Val!=0) g4_dyn += g4_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z4, 4);
      if (g4_prime6Val!=0) g4_dyn += g4_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z4, 4);
      if (g4_prime7Val!=0) g4_dyn += g4_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z4, 4);

      if (cz_q1sq!=0.) g4_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z41, 2));
      if (cz_q2sq!=0.) g4_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z42, 2));
      if (cz_q12sq!=0.) g4_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z40, 2));

      g1_dynIm = 0;
      if (g1_primeValIm!=0) g1_dynIm += g1_primeValIm * pow(Lambda_z1, 4)/(pow(Lambda_z1, 2) + pow(m1_, 2))/(pow(Lambda_z1, 2) + pow(m2_, 2));
      if (g1_prime2ValIm!=0) g1_dynIm += g1_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z1, 2);
      if (g1_prime3ValIm!=0) g1_dynIm += g1_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z1, 2);
      if (g1_prime4ValIm!=0) g1_dynIm += g1_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
      if (g1_prime5ValIm!=0) g1_dynIm += g1_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z1, 4);
      if (g1_prime6ValIm!=0) g1_dynIm += g1_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z1, 4);
      if (g1_prime7ValIm!=0) g1_dynIm += g1_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z1, 4);

      if (cz_q1sq!=0.) g1_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z11, 2));
      if (cz_q2sq!=0.) g1_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z12, 2));
      if (cz_q12sq!=0.) g1_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z10, 2));
      g1_dynIm += g1ValIm;

      g2_dynIm = g2ValIm;
      if (g2_primeValIm!=0) g2_dynIm += g2_primeValIm * pow(Lambda_z2, 4)/(pow(Lambda_z2, 2) + pow(m1_, 2))/(pow(Lambda_z2, 2) + pow(m2_, 2));
      if (g2_prime2ValIm!=0) g2_dynIm += g2_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z2, 2);
      if (g2_prime3ValIm!=0) g2_dynIm += g2_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z2, 2);
      if (g2_prime4ValIm!=0) g2_dynIm += g2_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
      if (g2_prime5ValIm!=0) g2_dynIm += g2_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z2, 4);
      if (g2_prime6ValIm!=0) g2_dynIm += g2_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z2, 4);
      if (g2_prime7ValIm!=0) g2_dynIm += g2_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z2, 4);

      if (cz_q1sq!=0.) g2_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z21, 2));
      if (cz_q2sq!=0.) g2_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z22, 2));
      if (cz_q12sq!=0.) g2_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z20, 2));

      g3_dynIm = g3ValIm;
      if (g3_primeValIm!=0) g3_dynIm += g3_primeValIm * pow(Lambda_z3, 4)/(pow(Lambda_z3, 2) + pow(m1_, 2))/(pow(Lambda_z3, 2) + pow(m2_, 2));
      if (g3_prime2ValIm!=0) g3_dynIm += g3_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z3, 2);
      if (g3_prime3ValIm!=0) g3_dynIm += g3_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z3, 2);
      if (g3_prime4ValIm!=0) g3_dynIm += g3_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
      if (g3_prime5ValIm!=0) g3_dynIm += g3_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z3, 4);
      if (g3_prime6ValIm!=0) g3_dynIm += g3_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z3, 4);
      if (g3_prime7ValIm!=0) g3_dynIm += g3_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z3, 4);

      if (cz_q1sq!=0.) g3_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z31, 2));
      if (cz_q2sq!=0.) g3_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z32, 2));
      if (cz_q12sq!=0.) g3_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z30, 2));

      g4_dynIm = g4ValIm;
      if (g4_primeValIm!=0) g4_dynIm += g4_primeValIm * pow(Lambda_z4, 4)/(pow(Lambda_z4, 2) + pow(m1_, 2))/(pow(Lambda_z4, 2) + pow(m2_, 2));
      if (g4_prime2ValIm!=0) g4_dynIm += g4_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z4, 2);
      if (g4_prime3ValIm!=0) g4_dynIm += g4_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z4, 2);
      if (g4_prime4ValIm!=0) g4_dynIm += g4_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
      if (g4_prime5ValIm!=0) g4_dynIm += g4_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z4, 4);
      if (g4_prime6ValIm!=0) g4_dynIm += g4_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z4, 4);
      if (g4_prime7ValIm!=0) g4_dynIm += g4_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z4, 4);

      if (cz_q1sq!=0.) g4_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z41, 2));
      if (cz_q2sq!=0.) g4_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z42, 2));
      if (cz_q12sq!=0.) g4_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z40, 2));
    }
    else if ((VGammaVpmode1==0 && VGammaVpmode2==2) || (VGammaVpmode1==2 && VGammaVpmode2==0)){ // VVp/VpV
      g1_dyn = gvvp1Val;
      g1_dynIm = gvvp1ValIm;
    }
    else if (VGammaVpmode1==2 && VGammaVpmode2==2){ // VpVp
      g1_dyn = gvpvp1Val;
      g1_dynIm = gvpvp1ValIm;
    }
  }
  else if (!(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell) && (VGammaVpmode1+VGammaVpmode2)==1){ // Zgs/Zg
    if (gzgs1_prime2Val!=0 && VGammaVpmode1==0) g1_dyn += gzgs1_prime2Val* pow(m1_, 2)/pow(Lambda_zgs1, 2);
    if (gzgs1_prime2Val!=0 && VGammaVpmode2==0) g1_dyn += gzgs1_prime2Val* pow(m2_, 2)/pow(Lambda_zgs1, 2);
    g2_dyn = gzgs2Val;
    g3_dyn = gzgs3Val;
    g4_dyn = gzgs4Val;

    if (gzgs1_prime2ValIm!=0 && VGammaVpmode1==0) g1_dynIm += gzgs1_prime2ValIm* pow(m1_, 2)/pow(Lambda_zgs1, 2);
    if (gzgs1_prime2ValIm!=0 && VGammaVpmode2==0) g1_dynIm += gzgs1_prime2ValIm* pow(m2_, 2)/pow(Lambda_zgs1, 2);
    g2_dynIm = gzgs2ValIm;
    g3_dynIm = gzgs3ValIm;
    g4_dynIm = gzgs4ValIm;
  }
  else{ // gsgs/gsg/gg
    g2_dyn = ggsgs2Val;
    g3_dyn = ggsgs3Val;
    g4_dyn = ggsgs4Val;

    g2_dynIm = ggsgs2ValIm;
    g3_dynIm = ggsgs3ValIm;
    g4_dynIm = ggsgs4ValIm;
  }
  /*
  MELAout << "g1: " << g1_dyn << " " << g1_dynIm << '\t';
  MELAout << "g2: " << g2_dyn << " " << g2_dynIm << '\t';
  MELAout << "g3: " << g3_dyn << " " << g3_dynIm << '\t';
  MELAout << "g4: " << g4_dyn << " " << g4_dynIm << '\t';
  MELAout << endl;
  */

  Double_t kappa = s/pow(Lambda, 2);
  a3Re = -2.*g4_dyn*pow(m12, 2);
  a2Re = -(2.*g2_dyn + g3_dyn*kappa)*pow(m12, 2);
  a1Re = g1_dyn*pow(mV, 2) - a2Re*s/pow(m12, 2);
  a3Im = -2.*g4_dynIm*pow(m12, 2);
  a2Im = -(2.*g2_dynIm + g3_dynIm*kappa)*pow(m12, 2);
  a1Im = g1_dynIm*pow(mV, 2) - a2Im*s/pow(m12, 2);

  if (verbosity>=TVar::DEBUG){
    MELAout << "RooSpinZero::calculateAi( "
      << a1Re << " , " << a1Im << " , " << a2Re << " , " << a2Im << " , " << a3Re << " , " << a3Im << " , " << VGammaVpmode1 << " , " << VGammaVpmode2
      << " ):\n"
      << "\t- Input:\n"
      << "\t- mV, m1, m2, m12, s = " << mV << " , " << m1 << " , " << m2 << " , " << m12 << " , " << s << '\n'
      << "\t- g1_dyn = (" << g1_dyn << " , " << g1_dynIm << ")" << '\n'
      << "\t- g2_dyn = (" << g2_dyn << " , " << g2_dynIm << ")" << '\n'
      << "\t- g3_dyn = (" << g3_dyn << " , " << g3_dynIm << ")" << '\n'
      << "\t- g4_dyn = (" << g4_dyn << " , " << g4_dynIm << ")"
      << endl;
  }
}
void RooSpinZero::calculateAmplitudes(
  Double_t& A00Re, Double_t& A00Im, Double_t& AppRe, Double_t& AppIm, Double_t& AmmRe, Double_t& AmmIm,
  int VGammaVpmode1, int VGammaVpmode2
)const{
  if (verbosity>=TVar::DEBUG) MELAout << "Begin RooSpinZero::calculateAmplitudes( " << VGammaVpmode1 << " , " << VGammaVpmode2 << " )" << endl;

  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;

  Double_t a1Re, a2Re, a3Re, a1Im, a2Im, a3Im;
  calculateAi(a1Re, a1Im, a2Re, a2Im, a3Re, a3Im, VGammaVpmode1, VGammaVpmode2);

  Double_t propV1Re=1, propV2Re=1, propHRe=1;
  Double_t propV1Im=0, propV2Im=0, propHIm=0;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell){
    int proptype=-1;
    switch (VGammaVpmode1){
    case 0:
      proptype=1;
      break;
    case 1:
      proptype=0;
      break;
    case 2:
      proptype=3;
      break;
    default:
      break;
    }
    calculatePropagator(propV1Re, propV1Im, m1_, proptype);
  }
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell){
    int proptype=-1;
    switch (VGammaVpmode2){
    case 0:
      proptype=1;
      break;
    case 1:
      proptype=0;
      break;
    case 2:
      proptype=3;
      break;
    default:
      break;
    }
    calculatePropagator(propV2Re, propV2Im, m2_, proptype);
  }
  calculatePropagator(propHRe, propHIm, m12, 2);

  Double_t ampScale = calculateAmplitudeScale(VGammaVpmode2, VGammaVpmode2)*pow(GeVunit, 2);

  Double_t eta1 = m1_ / m12;
  Double_t eta2 = m2_ / m12;
  Double_t eta1p2 = 1.;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) eta1p2 *= eta1;
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) eta1p2 *= eta2;

  // etas = s/m12**2
  Double_t etas = (1. - pow(eta1, 2) - pow(eta2, 2))/2.;
  if (pow(eta1+eta2, 2)>1.) etas = -etas;
  // x = (s/(m1*m2))**2 - 1 = (etas**2 - (eta1*eta2)**2)/(eta1*eta2)**2
  Double_t etasp = pow(etas, 2) - pow(eta1*eta2, 2); // Notice how eta1p2 is not used. The second set of multiplications below is the reason to it!
  if (etasp<0) etasp=0;

  if (verbosity>=TVar::DEBUG) MELAout << "RooSpinZero::calculateAmplitudes:\n"
    << "\t- eta1, eta2, eta1p2, etasp = "
    << eta1 << " , " << eta2 << " , " << eta1p2 << " , " << etasp << '\n'
    << "\t- ampScale, prop1, prop2 = "
    << ampScale << " , " << propV1Re << "," << propV1Im << " , " << propV2Re << "," << propV2Im
    << endl;

  Double_t A00Re_tmp=0, A00Im_tmp=0, AppRe_tmp=0, AppIm_tmp=0, AmmRe_tmp=0, AmmIm_tmp=0;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell){
    A00Re_tmp = -2.*(a1Re*etas + a2Re*etasp)*ampScale;
    A00Im_tmp = -2.*(a1Im*etas + a2Im*etasp)*ampScale;
  }
  AppRe_tmp = (a1Re - a3Im*sqrt(etasp))*ampScale;
  AppIm_tmp = (a1Im + a3Re*sqrt(etasp))*ampScale;
  AmmRe_tmp = (a1Re + a3Im*sqrt(etasp))*ampScale;
  AmmIm_tmp = (a1Im - a3Re*sqrt(etasp))*ampScale;

  AppRe_tmp *= eta1p2;
  AppIm_tmp *= eta1p2;
  AmmRe_tmp *= eta1p2;
  AmmIm_tmp *= eta1p2;

  if (verbosity>=TVar::DEBUG) MELAout << "RooSpinZero::calculateAmplitudes:\n"
    << "\t- A00 factor = " << A00Re_tmp << " , " << A00Im_tmp << '\n'
    << "\t- A++ factor = " << AppRe_tmp << " , " << AppIm_tmp << '\n'
    << "\t- A-- factor = " << AmmRe_tmp << " , " << AmmIm_tmp
    << endl;

  // A_old = ARe_old+i*AIm_old => A_new = ARe_new + i*AIm_new = A_old*propV1*propV2
  std::vector<Double_t> A00_reals, A00_imags; A00_reals.push_back(A00Re_tmp); A00_imags.push_back(A00Im_tmp); A00_reals.push_back(propV1Re); A00_imags.push_back(propV1Im); A00_reals.push_back(propV2Re); A00_imags.push_back(propV2Im); A00_reals.push_back(propHRe); A00_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(A00_reals, A00_imags, A00Re, A00Im);
  std::vector<Double_t> App_reals, App_imags; App_reals.push_back(AppRe_tmp); App_imags.push_back(AppIm_tmp); App_reals.push_back(propV1Re); App_imags.push_back(propV1Im); App_reals.push_back(propV2Re); App_imags.push_back(propV2Im); App_reals.push_back(propHRe); App_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(App_reals, App_imags, AppRe, AppIm);
  std::vector<Double_t> Amm_reals, Amm_imags; Amm_reals.push_back(AmmRe_tmp); Amm_imags.push_back(AmmIm_tmp); Amm_reals.push_back(propV1Re); Amm_imags.push_back(propV1Im); Amm_reals.push_back(propV2Re); Amm_imags.push_back(propV2Im); Amm_reals.push_back(propHRe); Amm_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(Amm_reals, Amm_imags, AmmRe, AmmIm);

  if (!(A00Re==A00Re) || !(AppRe==AppRe) || !(AmmRe==AmmRe) || !(A00Im==A00Im) || !(AppIm==AppIm) || !(AmmIm==AmmIm)){
    MELAout << "A00 = " << A00Re << ", " << A00Im << endl;
    MELAout << "A++ = " << AppRe << ", " << AppIm << endl;
    MELAout << "A-- = " << AmmRe << ", " << AmmIm << endl;
    MELAout << eta1 << '\t' << eta2 << '\t' << etas << '\t' << etasp << '\t' << eta1p2 << endl;
  }
}

Bool_t RooSpinZero::computeNeededAmplitude(int VGammaVpmode1, int VGammaVpmode2) const{
  return (
    (
    (VGammaVpmode1==VGammaVpmode2 && VGammaVpmode1==0) && (
      g1Val!=0.
      || g1_primeVal!=0.
      || g1_prime2Val!=0.
      || g1_prime3Val!=0.
      || g1_prime4Val!=0.
      || g1_prime5Val!=0.
      || g1_prime6Val!=0.
      || g1_prime7Val!=0.
      || g2Val!=0.
      || g2_primeVal!=0.
      || g2_prime2Val!=0.
      || g2_prime3Val!=0.
      || g2_prime4Val!=0.
      || g2_prime5Val!=0.
      || g2_prime6Val!=0.
      || g2_prime7Val!=0.
      || g4Val!=0.
      || g4_primeVal!=0.
      || g4_prime2Val!=0.
      || g4_prime3Val!=0.
      || g4_prime4Val!=0.
      || g4_prime5Val!=0.
      || g4_prime6Val!=0.
      || g4_prime7Val!=0.
      || g3Val!=0.
      || g3_primeVal!=0.
      || g3_prime2Val!=0.
      || g3_prime3Val!=0.
      || g3_prime4Val!=0.
      || g3_prime5Val!=0.
      || g3_prime6Val!=0.
      || g3_prime7Val!=0.

      || g1ValIm!=0.
      || g1_primeValIm!=0.
      || g1_prime2ValIm!=0.
      || g1_prime3ValIm!=0.
      || g1_prime4ValIm!=0.
      || g1_prime5ValIm!=0.
      || g1_prime6ValIm!=0.
      || g1_prime7ValIm!=0.
      || g2ValIm!=0.
      || g2_primeValIm!=0.
      || g2_prime2ValIm!=0.
      || g2_prime3ValIm!=0.
      || g2_prime4ValIm!=0.
      || g2_prime5ValIm!=0.
      || g2_prime6ValIm!=0.
      || g2_prime7ValIm!=0.
      || g4ValIm!=0.
      || g4_primeValIm!=0.
      || g4_prime2ValIm!=0.
      || g4_prime3ValIm!=0.
      || g4_prime4ValIm!=0.
      || g4_prime5ValIm!=0.
      || g4_prime6ValIm!=0.
      || g4_prime7ValIm!=0.
      || g3ValIm!=0.
      || g3_primeValIm!=0.
      || g3_prime2ValIm!=0.
      || g3_prime3ValIm!=0.
      || g3_prime4ValIm!=0.
      || g3_prime5ValIm!=0.
      || g3_prime6ValIm!=0.
      || g3_prime7ValIm!=0.
      )
      )
    ||
    (
    ((VGammaVpmode1==0 && VGammaVpmode2==1) || (VGammaVpmode2==0 && VGammaVpmode1==1)) && (
      gzgs1_prime2Val!=0.
      || gzgs2Val!=0.
      || gzgs4Val!=0.
      || gzgs3Val!=0.

      || gzgs1_prime2ValIm!=0.
      || gzgs2ValIm!=0.
      || gzgs4ValIm!=0.
      || gzgs3ValIm!=0.
      )
      )
    ||
    (
    (VGammaVpmode1==VGammaVpmode2 && VGammaVpmode1==1) && (
      ggsgs2Val!=0.
      || ggsgs4Val!=0.
      || ggsgs3Val!=0.

      || ggsgs2ValIm!=0.
      || ggsgs4ValIm!=0.
      || ggsgs3ValIm!=0.
      )
      )
    ||
    (
    ((VGammaVpmode1==0 && VGammaVpmode2==2) || (VGammaVpmode2==0 && VGammaVpmode1==2)) && (
      gvvp1Val!=0.

      || gvvp1ValIm!=0.
      )
      )
    ||
    (
    (VGammaVpmode1==VGammaVpmode2 && VGammaVpmode1==2) && (
      gvpvp1Val!=0.

      || gvpvp1ValIm!=0.
      )
      )
    );
}

void RooSpinZero::printParameters() const{
  MELAout << "g1: (" << g1Val << ", " << g1ValIm << ")" << endl;
  MELAout << "g2: (" << g2Val << ", " << g2ValIm << ")" << endl;
  MELAout << "g3: (" << g3Val << ", " << g3ValIm << ")" << endl;
  MELAout << "g4: (" << g4Val << ", " << g4ValIm << ")" << endl;
  MELAout << "g1_prime: (" << g1_primeVal << ", " << g1_primeValIm << ")" << endl;
  MELAout << "g2_prime: (" << g2_primeVal << ", " << g2_primeValIm << ")" << endl;
  MELAout << "g3_prime: (" << g3_primeVal << ", " << g3_primeValIm << ")" << endl;
  MELAout << "g4_prime: (" << g4_primeVal << ", " << g4_primeValIm << ")" << endl;
  MELAout << "g1_prime2: (" << g1_prime2Val << ", " << g1_prime2ValIm << ")" << endl;
  MELAout << "g2_prime2: (" << g2_prime2Val << ", " << g2_prime2ValIm << ")" << endl;
  MELAout << "g3_prime2: (" << g3_prime2Val << ", " << g3_prime2ValIm << ")" << endl;
  MELAout << "g4_prime2: (" << g4_prime2Val << ", " << g4_prime2ValIm << ")" << endl;
  MELAout << "g1_prime3: (" << g1_prime3Val << ", " << g1_prime3ValIm << ")" << endl;
  MELAout << "g2_prime3: (" << g2_prime3Val << ", " << g2_prime3ValIm << ")" << endl;
  MELAout << "g3_prime3: (" << g3_prime3Val << ", " << g3_prime3ValIm << ")" << endl;
  MELAout << "g4_prime3: (" << g4_prime3Val << ", " << g4_prime3ValIm << ")" << endl;
  MELAout << "g1_prime4: (" << g1_prime4Val << ", " << g1_prime4ValIm << ")" << endl;
  MELAout << "g2_prime4: (" << g2_prime4Val << ", " << g2_prime4ValIm << ")" << endl;
  MELAout << "g3_prime4: (" << g3_prime4Val << ", " << g3_prime4ValIm << ")" << endl;
  MELAout << "g4_prime4: (" << g4_prime4Val << ", " << g4_prime4ValIm << ")" << endl;
  MELAout << "g1_prime5: (" << g1_prime5Val << ", " << g1_prime5ValIm << ")" << endl;
  MELAout << "g2_prime5: (" << g2_prime5Val << ", " << g2_prime5ValIm << ")" << endl;
  MELAout << "g3_prime5: (" << g3_prime5Val << ", " << g3_prime5ValIm << ")" << endl;
  MELAout << "g4_prime5: (" << g4_prime5Val << ", " << g4_prime5ValIm << ")" << endl;
  MELAout << "g1_prime6: (" << g1_prime6Val << ", " << g1_prime6ValIm << ")" << endl;
  MELAout << "g2_prime6: (" << g2_prime6Val << ", " << g2_prime6ValIm << ")" << endl;
  MELAout << "g3_prime6: (" << g3_prime6Val << ", " << g3_prime6ValIm << ")" << endl;
  MELAout << "g4_prime6: (" << g4_prime6Val << ", " << g4_prime6ValIm << ")" << endl;
  MELAout << "g1_prime7: (" << g1_prime7Val << ", " << g1_prime7ValIm << ")" << endl;
  MELAout << "g2_prime7: (" << g2_prime7Val << ", " << g2_prime7ValIm << ")" << endl;
  MELAout << "g3_prime7: (" << g3_prime7Val << ", " << g3_prime7ValIm << ")" << endl;
  MELAout << "g4_prime7: (" << g4_prime7Val << ", " << g4_prime7ValIm << ")" << endl;
  MELAout << "gzgs1_prime2: (" << gzgs1_prime2Val << ", " << gzgs1_prime2ValIm << ")" << endl;
  MELAout << "gzgs2: (" << gzgs2Val << ", " << gzgs2ValIm << ")" << endl;
  MELAout << "gzgs3: (" << gzgs3Val << ", " << gzgs3ValIm << ")" << endl;
  MELAout << "gzgs4: (" << gzgs4Val << ", " << gzgs4ValIm << ")" << endl;
  MELAout << "ggsgs2: (" << ggsgs2Val << ", " << ggsgs2ValIm << ")" << endl;
  MELAout << "ggsgs3: (" << ggsgs3Val << ", " << ggsgs3ValIm << ")" << endl;
  MELAout << "ggsgs4: (" << ggsgs4Val << ", " << ggsgs4ValIm << ")" << endl;

  MELAout << "Lambda: " << Lambda << endl;
  MELAout << "Lambda_zgs1: " << Lambda_zgs1 << endl;
  MELAout << "Lambda_z1: " << Lambda_z1 << endl;
  MELAout << "Lambda_z2: " << Lambda_z2 << endl;
  MELAout << "Lambda_z3: " << Lambda_z3 << endl;
  MELAout << "Lambda_z4: " << Lambda_z4 << endl;
  MELAout << "Lambda_Q: " << Lambda_Q << endl;
  MELAout << "Lambda_z11: " << Lambda_z11 << endl;
  MELAout << "Lambda_z21: " << Lambda_z21 << endl;
  MELAout << "Lambda_z31: " << Lambda_z31 << endl;
  MELAout << "Lambda_z41: " << Lambda_z41 << endl;
  MELAout << "Lambda_z12: " << Lambda_z12 << endl;
  MELAout << "Lambda_z22: " << Lambda_z22 << endl;
  MELAout << "Lambda_z32: " << Lambda_z32 << endl;
  MELAout << "Lambda_z42: " << Lambda_z42 << endl;
  MELAout << "Lambda_z10: " << Lambda_z10 << endl;
  MELAout << "Lambda_z20: " << Lambda_z20 << endl;
  MELAout << "Lambda_z30: " << Lambda_z30 << endl;
  MELAout << "Lambda_z40: " << Lambda_z40 << endl;
  MELAout << "cz_q1sq: " << cz_q1sq << endl;
  MELAout << "cz_q2sq: " << cz_q2sq << endl;
  MELAout << "cz_q12sq: " << cz_q12sq << endl;

  MELAout << "gvvp1: (" << gvvp1Val << ", " << gvvp1ValIm << ")" << endl;
  MELAout << "gvpvp1: (" << gvpvp1Val << ", " << gvpvp1ValIm << ")" << endl;

  RooSpin::printParameters();
}
