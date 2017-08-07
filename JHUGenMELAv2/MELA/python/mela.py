"""
Python wrapper for MELA
>>> from mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar
>>> m = Mela(13, 125)
Then you can use m mostly like a C++ mela object.

>>> m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)

The main change is that computeP and similar functions return the ME instead of modifying a reference.
(computeP_selfD* are not implemented here, but you can access the functionality through computeP)

>>> print m.computeP(False)

You can modify the couplings either by
>>> m.selfDHzzcoupl[0][0][0] = 1
>>> m.selfDHzzcoupl[0][0][1] = 2
or by the more convenient way
>>> m.ghz1 = 1+2j

SimpleParticle_t and SimpleParticleCollection_t are implemented to take inputs in ways more normal for python.
>>> daughters = SimpleParticleCollection_t([
...                                         "11 -71.89 30.50 -47.20 91.25",
...                                         #...other daughters
...                                        ])

There's also a function for convenience
>>> m.setInputEvent_fromLHE('''
... <event>
... #(...)
... </event>
... ''')
which is useful for quick tests.

See examples at the bottom.
"""

from collections import namedtuple
import ROOT
from pythonmelautils import MultiDimensionalCppArray, NamedTemporaryMacro, SelfDParameter, SelfDCoupling
from ROOT import TUtil, TVar

class Mela(object):
  counter = 0
  doneinit = False
  computeptemplate = """
    #include <Mela.h>
    float getPAux(Mela& mela) {
      float result;
      mela.getPAux(result);
      return result;
    }
    vector<float> computeDecayAngles(Mela& mela) {
      vector<float> result(8);
      mela.computeDecayAngles(
        result[0],
        result[1],
        result[2],
        result[3],
        result[4],
        result[5],
        result[6],
        result[7]
      );
      return result;
    }
    //not implementing the computeP_selfD* functions here
    //would be easier to do in pure python but not worth it anyway
    float computeP(Mela& mela, bool useConstant) {
      float result;
      mela.computeP(result, useConstant);
      return result;
    }
    float computeD_CP(Mela& mela, TVar::MatrixElement myME, TVar::Process myType) {
      float result;
      mela.computeD_CP(myME, myType, result);
      return result;
    }
    float computeProdP(Mela& mela, bool useConstant) {
      float result;
      mela.computeProdP(result, useConstant);
      return result;
    }
    float computeProdDecP(Mela& mela, bool useConstant) {
      float result;
      mela.computeProdDecP(result, useConstant);
      return result;
    }
    float computeProdP_VH(Mela& mela, bool includeHiggsDecay, bool useConstant) {
      float result;
      mela.computeProdP_VH(result, includeHiggsDecay, useConstant);
      return result;
    }
    float computeProdP_ttH(Mela& mela, int topProcess, int topDecay, bool useConstant) {
      float result;
      mela.computeProdP_ttH(result, topProcess, topDecay, useConstant);
      return result;
    }
    float compute4FermionWeight(Mela& mela) {
      float result;
      mela.compute4FermionWeight(result);
      return result;
    }
    float getXPropagator(Mela& mela, TVar::ResonancePropagatorScheme scheme) {
      float result;
      mela.getXPropagator(scheme, result);
      return result;
    }
    float computePM4l(Mela& mela, TVar::SuperMelaSyst syst) {
      float result;
      mela.computePM4l(syst, result);
      return result;
    }
    float computeD_gg(Mela& mela, TVar::MatrixElement myME, TVar::Process myType) {
      float result;
      mela.computeD_gg(myME, myType, result);
      return result;
    }
    float getConstant(Mela& mela) {
      float result;
      mela.getConstant(result);
      return result;
    }
    float computeDijetConvBW(Mela& mela) {
      float result;
      mela.computeDijetConvBW(result);
      return result;
    }
  """
  def __init__(self, *args, **kwargs):
    self.__mela = ROOT.Mela(*args, **kwargs)
    self.index = self.counter
    type(self).counter += 1

    arrays  = (
               ("selfDHggcoupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HGG, 2)),
               ("selfDHg4g4coupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HGG, 2)),
               ("selfDHqqcoupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HQQ, 2)),
               ("selfDHbbcoupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HQQ, 2)),
               ("selfDHttcoupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HQQ, 2)),
               ("selfDHb4b4coupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HQQ, 2)),
               ("selfDHt4t4coupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HQQ, 2)),
               ("selfDHzzcoupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HVV, 2)),
               ("selfDHwwcoupl", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HVV, 2)),
               ("selfDHzzLambda_qsq", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HVV_LAMBDAQSQ, ROOT.py_SIZE_HVV_CQSQ)),
               ("selfDHwwLambda_qsq", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HVV_LAMBDAQSQ, ROOT.py_SIZE_HVV_CQSQ)),
               ("selfDHzzCLambda_qsq", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HVV_CQSQ)),
               ("selfDHwwCLambda_qsq", (ROOT.nSupportedHiggses, ROOT.py_SIZE_HVV_CQSQ)),
               ("selfDHzzpcoupl", (ROOT.py_SIZE_HVV, 2)),
               ("selfDHzpzpcoupl", (ROOT.py_SIZE_HVV, 2)),
               ("selfDHzpcontact", (ROOT.py_SIZE_Vp, 2)),
               ("selfDHwwpcoupl", (ROOT.py_SIZE_HVV, 2)),
               ("selfDHwpwpcoupl", (ROOT.py_SIZE_HVV, 2)),
               ("selfDHwpcontact", (ROOT.py_SIZE_Vp, 2)),
               ("selfDZqqcoupl", (ROOT.py_SIZE_ZQQ, 2)),
               ("selfDZvvcoupl", (ROOT.py_SIZE_ZVV, 2)),
               ("selfDGqqcoupl", (ROOT.py_SIZE_GQQ, 2)),
               ("selfDGggcoupl", (ROOT.py_SIZE_GGG, 2)),
               ("selfDGvvcoupl", (ROOT.py_SIZE_GVV, 2)),
              )

    f = None
    for name, dimensions in arrays:
        setattr(
                self,
                name,
                MultiDimensionalCppArray(
                                         "mela{}{}".format(self.index, name),
                                         "mela.{}".format(name),
                                         ["Mela.h"],
                                         {"Mela& mela": self.__mela},
                                         *dimensions
                                        )
               )
        f = getattr(self, name).writecpp(f)
    bkpgErrorIgnoreLevel, ROOT.gErrorIgnoreLevel = ROOT.gErrorIgnoreLevel, ROOT.kInfo+1
    f.write(self.computeptemplate)
    for name, dimensions in arrays:
        getattr(self, name).compilecpp(f)
    ROOT.gErrorIgnoreLevel = bkpgErrorIgnoreLevel
    self.doneinit = True

  def __getattr__(self, name):
    return getattr(self.__mela, name)

  def __setattr__(self, name, value):
    if self.doneinit:
      if hasattr(self.__mela, name):
        return setattr(self.__mela, name, value)
      elif hasattr(type(self), name) and hasattr(getattr(type(self), name), "__set__"):
        return getattr(type(self), name).__set__(self, value)
      else:
        raise ValueError("Can't set attribute '{}' for python Mela, which doesn't exist for C++ MELA".format(name))
    else:
      super(Mela, self).__setattr__(name, value)

  def setInputEvent_fromLHE(self, event, isgen):
    lines = event.split("\n")
    lines = [line for line in lines if not ("<event>" in line or "</event>" in line or not line.split("#")[0].strip())]
    nparticles, _, _, _, _, _ = lines[0].split()
    nparticles = int(nparticles)
    if nparticles != len(lines)-1:
      raise ValueError("Wrong number of particles! Should be {}, have {}".replace(nparticles, len(lines)-1))
    daughters, mothers, associated = [], [], []
    ids = [None]
    mother1s = [None]
    mother2s = [None]
    for line in lines[1:]:
      id, status, mother1, mother2 = (int(_) for _ in line.split()[0:4])
      ids.append(id)
      mother1s.append(mother1)
      mother2s.append(mother2)
      if (1 <= abs(id) <= 6 or abs(id) == 21) and not isgen:
        line = line.replace(str(id), "0", 1)  #replace the first instance of the jet id with 0, which means unknown jet
      if status == -1:
        mothers.append(line)
      elif status == 1 and (1 <= abs(id) <= 6 or 11 <= abs(id) <= 16 or abs(id) in (21, 22)):
        while True:
          if mother1 != mother2 or mother1 is None:
            associated.append(line)
            break
          if ids[mother1] == 25:
            daughters.append(line)
            break
          mother2 = mother2s[mother1]
          mother1 = mother1s[mother1]
    #print "mothers"
    #for _ in mothers: print _
    #print "daughters"
    #for _ in daughters: print _
    #print "associated"
    #for _ in associated: print _
    self.setInputEvent(SimpleParticleCollection_t(daughters), SimpleParticleCollection_t(associated), SimpleParticleCollection_t(mothers), isgen)

  def getPAux(self): return ROOT.getPAux(self.__mela)
  DecayAngles = namedtuple("DecayAngles", "qH m1 m2 costheta1 costheta2 Phi costhetastar Phi1")
  def computeDecayAngles(self): return self.DecayAngles(*ROOT.computeDecayAngles(self.__mela))
  def computeP(self, useConstant=True): return ROOT.computeP(self.__mela, useConstant)
  def computeD_CP(self, myME, myType): return ROOT.computeD_CP(self.__mela, myME, myType)
  def computeProdP(self, useConstant=True): return ROOT.computeProdP(self.__mela, useConstant)
  def computeProdDecP(self, useConstant=True): return ROOT.computeProdDecP(self.__mela, useConstant)
  def compute4FermionWeight(self): return ROOT.compute4FermionWeight(self.__mela)
  def getXPropagator(self, scheme): return ROOT.getXPropagator(self.__mela, scheme)
  def computePM4l(self, syst): return ROOT.computePM4l(self.__mela, syst)
  def computeD_gg(self, myME, myType): return ROOT.computeD_gg(self.__mela, myME, myType)
  def computeProdP_VH(self, includeHiggsDecay=False, useConstant=True): return ROOT.computeProdP_VH(self.__mela, includeHiggsDecay, useConstant)
  def computeProdP_ttH(self, topProcess=2, topDecay=0, useConstant=True): return ROOT.computeProdP_ttH(self.__mela, topProcess, topDecay, useConstant)
  def getConstant(self): return ROOT.getConstant(self.__mela)
  def computeDijetConvBW(self): return ROOT.computeDijetConvBW(self.__mela)

  ghg2 = SelfDCoupling("selfDHggcoupl", 0, ROOT.py_gHIGGS_GG_2)
  ghg3 = SelfDCoupling("selfDHggcoupl", 0, ROOT.py_gHIGGS_GG_3)
  ghg4 = SelfDCoupling("selfDHggcoupl", 0, ROOT.py_gHIGGS_GG_4)

  #https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/41232f911b4f03065ae2b83752b5bcd4daacaa2c/MELA/fortran/mod_JHUGenMELA.F90#L123-L168
  ghz1 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1)
  ghz2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2)
  ghz3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3)
  ghz4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_4)

  ghzgs2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_ZA_2)
  ghzgs3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_ZA_3)
  ghzgs4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_ZA_4)
  ghgsgs2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_AA_2)
  ghgsgs3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_AA_3)
  ghgsgs4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_AA_4)

  ghz1_prime = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME)
  ghz1_prime2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME2)
  ghz1_prime3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME3)
  ghz1_prime4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME4)
  ghz1_prime5 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME5)

  ghz2_prime = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME)
  ghz2_prime2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME2)
  ghz2_prime3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME3)
  ghz2_prime4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME4)
  ghz2_prime5 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME5)

  ghz3_prime = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME)
  ghz3_prime2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME2)
  ghz3_prime3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME3)
  ghz3_prime4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME4)
  ghz3_prime5 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME5)

  ghz4_prime = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME)
  ghz4_prime2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME2)
  ghz4_prime3 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME3)
  ghz4_prime4 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME4)
  ghz4_prime5 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME5)

  ghzgs1_prime2 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_ZA_1_PRIME2)

  ghz1_prime6 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME6)
  ghz1_prime7 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME7)
  ghz2_prime6 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME6)
  ghz2_prime7 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME7)
  ghz3_prime6 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME6)
  ghz3_prime7 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME7)
  ghz4_prime6 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME6)
  ghz4_prime7 = SelfDCoupling("selfDHzzcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME7)

  cz_q1sq = SelfDParameter("selfDHzzCLambda_qsq", 0, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_z11 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_1, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_z12 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_2, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_z13 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_3, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_z14 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_4, ROOT.py_cLambdaHIGGS_VV_QSQ1)

  cz_q2sq = SelfDParameter("selfDHzzCLambda_qsq", 0, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_z21 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_1, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_z22 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_2, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_z23 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_3, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_z24 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_4, ROOT.py_cLambdaHIGGS_VV_QSQ2)

  cz_q12sq = SelfDParameter("selfDHzzCLambda_qsq", 0, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_z01 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_1, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_z02 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_2, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_z03 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_3, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_z04 = SelfDParameter("selfDHzzLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_4, ROOT.py_cLambdaHIGGS_VV_QSQ12)

  ghw1 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1)
  ghw2 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2)
  ghw3 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3)
  ghw4 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_4)

  ghw1_prime = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME)
  ghw1_prime2 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME2)
  ghw1_prime3 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME3)
  ghw1_prime4 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME4)
  ghw1_prime5 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME5)

  ghw2_prime = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME)
  ghw2_prime2 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME2)
  ghw2_prime3 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME3)
  ghw2_prime4 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME4)
  ghw2_prime5 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME5)

  ghw3_prime = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME)
  ghw3_prime2 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME2)
  ghw3_prime3 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME3)
  ghw3_prime4 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME4)
  ghw3_prime5 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME5)

  ghw4_prime = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME)
  ghw4_prime2 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME2)
  ghw4_prime3 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME3)
  ghw4_prime4 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME4)
  ghw4_prime5 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_4_PRIME5)

  ghw1_prime6 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME6)
  ghw1_prime7 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_1_PRIME7)
  ghw2_prime6 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME6)
  ghw2_prime7 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_2_PRIME7)
  ghw3_prime6 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME6)
  ghw3_prime7 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME7)
  ghw4_prime6 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME6)
  ghw4_prime7 = SelfDCoupling("selfDHwwcoupl", 0, ROOT.py_gHIGGS_VV_3_PRIME7)

  cw_q1sq = SelfDParameter("selfDHwwCLambda_qsq", 0, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_w11 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_1, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_w12 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_2, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_w13 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_3, ROOT.py_cLambdaHIGGS_VV_QSQ1)
  Lambda_w14 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_4, ROOT.py_cLambdaHIGGS_VV_QSQ1)

  cw_q2sq = SelfDParameter("selfDHwwCLambda_qsq", 0, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_w21 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_1, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_w22 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_2, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_w23 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_3, ROOT.py_cLambdaHIGGS_VV_QSQ2)
  Lambda_w24 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_4, ROOT.py_cLambdaHIGGS_VV_QSQ2)

  cw_q12sq = SelfDParameter("selfDHwwCLambda_qsq", 0, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_w01 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_1, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_w02 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_2, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_w03 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_3, ROOT.py_cLambdaHIGGS_VV_QSQ12)
  Lambda_w04 = SelfDParameter("selfDHwwLambda_qsq", 0, ROOT.py_LambdaHIGGS_QSQ_VV_4, ROOT.py_cLambdaHIGGS_VV_QSQ12)

  kappa = SelfDCoupling("selfDHqqcoupl", 0, ROOT.py_gHIGGS_KAPPA)
  kappa_tilde = SelfDCoupling("selfDHqqcoupl", 0, ROOT.py_gHIGGS_KAPPA_TILDE)

  ghzzp1 = SelfDCoupling("selfDHzzpcoupl", ROOT.py_gHIGGS_VV_1)
  ghzpzp1 = SelfDCoupling("selfDHzpzpcoupl", ROOT.py_gHIGGS_VV_1)
  ezp_L_E = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_E)
  ezp_R_E = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_E)
  ezp_L_M = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_M)
  ezp_R_M = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_M)
  ezp_L_T = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_T)
  ezp_R_T = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_T)
  ezp_L_N = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_N)
  ezp_R_N = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_N)
  ezp_L_U = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_U)
  ezp_R_U = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_U)
  ezp_L_D = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_D)
  ezp_R_D = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_D)
  ezp_L_S = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_S)
  ezp_R_S = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_S)
  ezp_L_C = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_C)
  ezp_R_C = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_C)
  ezp_L_B = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_L_B)
  ezp_R_B = SelfDCoupling("selfDHzpcontact", ROOT.py_gHIGGS_Vp_R_B)

  ghwwp1 = SelfDCoupling("selfDHwwpcoupl", ROOT.py_gHIGGS_VV_1)
  ghwpwp1 = SelfDCoupling("selfDHwpwpcoupl", ROOT.py_gHIGGS_VV_1)
  ewp_L_E = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_L_E)
  ewp_R_E = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_R_E)
  ewp_L_M = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_L_M)
  ewp_R_M = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_R_M)
  ewp_L_T = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_L_T)
  ewp_R_T = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_R_T)
  ewp_L_U = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_L_U)
  ewp_R_U = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_R_U)
  ewp_L_C = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_L_C)
  ewp_R_C = SelfDCoupling("selfDHwpcontact", ROOT.py_gHIGGS_Vp_R_C)

  UseVprime = SelfDParameter("selfDUseVprime")
  M_Vprime = SelfDParameter("selfDM_Vprime")
  Ga_Vprime = SelfDParameter("selfDGa_Vprime")

  zprime_qq_left = SelfDCoupling("selfDZqqcoupl", ROOT.py_gZPRIME_QQ_LEFT)
  zprime_qq_right = SelfDCoupling("selfDZqqcoupl", ROOT.py_gZPRIME_QQ_RIGHT)
  zprime_zz_1 = SelfDCoupling("selfDZvvcoupl", ROOT.py_gZPRIME_VV_1)
  zprime_zz_2 = SelfDCoupling("selfDZvvcoupl", ROOT.py_gZPRIME_VV_2)

  graviton_qq_left = SelfDCoupling("selfDGqqcoupl", ROOT.py_gGRAVITON_QQ_LEFT)
  graviton_qq_right = SelfDCoupling("selfDGqqcoupl", ROOT.py_gGRAVITON_QQ_RIGHT)

  a1 = SelfDCoupling("selfDGggcoupl", ROOT.py_gGRAVITON_GG_1)
  a2 = SelfDCoupling("selfDGggcoupl", ROOT.py_gGRAVITON_GG_2)
  a3 = SelfDCoupling("selfDGggcoupl", ROOT.py_gGRAVITON_GG_3)
  a4 = SelfDCoupling("selfDGggcoupl", ROOT.py_gGRAVITON_GG_4)
  a5 = SelfDCoupling("selfDGggcoupl", ROOT.py_gGRAVITON_GG_5)

  b1 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_1)
  b2 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_2)
  b3 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_3)
  b4 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_4)
  b5 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_5)
  b6 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_6)
  b7 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_7)
  b8 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_8)
  b9 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_9)
  b10 = SelfDCoupling("selfDGvvcoupl", ROOT.py_gGRAVITON_VV_10)



def SimpleParticleCollection_t(iterable=None):
  if iterable is None: return ROOT.SimpleParticleCollection_t()
  result = ROOT.SimpleParticleCollection_t()
  for _ in iterable:
    result.push_back(SimpleParticle_t(_))
  return result

def SimpleParticle_t(lineorid, pxortlv=None, py=None, pz=None, e=None):
  if pxortlv is py is pz is e is None:
    if isinstance(lineorid, ROOT.SimpleParticle_t): return lineorid
    if isinstance(lineorid, basestring):
      lineorid = lineorid.split()
    if len(lineorid) == 13:
      id, status, mother1, mother2, color1, color2, px, py, pz, e, m, lifetime, spin = (f(_) for f, _ in zip((int, int, int, int, int, int, float, float, float, float, float, float, float), lineorid))
    elif len(lineorid) == 5:
      id, px, py, pz, e = (f(_) for f, _ in zip((int, float, float, float, float), lineorid))
    elif len(lineorid) == 2:
      pxortlv = lineorid[1]
    else:
      raise ValueError("len(lineorid) has to be 5 or 13, not {}".format(len(lineorid)))
  else:
    id = lineorid
    px = pxortlv

  if py is pz is e is None:
    tlv = pxortlv
  else:
    tlv = ROOT.TLorentzVector(px, py, pz, e)

  return ROOT.SimpleParticle_t(id, tlv)

if __name__ == "__main__":
  m = Mela()
  event1 = """
<event>
12  50   1.0000000E+00   1.2500000E+02   7.8125000E-03   1.2380607E-01
        2   -1    0    0  503    0  0.00000000000E+00  0.00000000000E+00  1.65430825479E+03  1.65430825479E+03  0.00000000000E+00 0.00000000000E+00  1.
       -1   -1    0    0    0  503  0.00000000000E+00  0.00000000000E+00 -1.42857195999E+01  1.42857195999E+01  0.00000000000E+00 0.00000000000E+00  1.
       24    2    1    2    0    0 -2.29473352103E+00 -1.04804828491E+02  4.95291431378E+02  5.12296652853E+02  7.83896718395E+01 0.00000000000E+00  1.
       25    2    1    2    0    0  2.29473352103E+00  1.04804828491E+02  1.14473110382E+03  1.15629732154E+03  1.24999511524E+02 0.00000000000E+00  1.
       14    1    3    3    0    0  4.42035961901E+00 -5.60456350211E+01  4.09886160671E+02  4.13723721213E+02  8.42936970218E-06 0.00000000000E+00  1.
      -13    1    3    3    0    0 -6.71509314004E+00 -4.87591934698E+01  8.54052707068E+01  9.85729316407E+01  1.05660000144E-01 0.00000000000E+00  1.
       23    2    4    4    0    0 -2.00748771644E+01  3.21702667586E+01  3.27018956548E+02  3.30034988785E+02  2.33188576920E+01 0.00000000000E+00  1.
       23    2    4    4    0    0  2.23696106855E+01  7.26345617324E+01  8.17712147272E+02  8.26262332755E+02  9.09950970840E+01 0.00000000000E+00  1.
      -11    1    7    7    0    0 -1.74223737299E+01  9.11950220870E+00  1.06644211152E+02  1.08442114510E+02  5.11001208360E-04 0.00000000000E+00  1.
       11    1    7    7    0    0 -2.65250343458E+00  2.30507645499E+01  2.20374745396E+02  2.21592874275E+02  5.10994690391E-04 0.00000000000E+00  1.
      -13    1    8    8    0    0  8.81223774828E+00  8.87930337607E+01  5.03683096793E+02  5.11525690007E+02  1.05660000328E-01 0.00000000000E+00  1.
       13    1    8    8    0    0  1.35573729372E+01 -1.61584720283E+01  3.14029050479E+02  3.14736642748E+02  1.05659999907E-01 0.00000000000E+00  1.
</event>
  """
  event2 = """
<event>
12  50   1.0000000E+00   1.2500000E+02   7.8125000E-03   1.2380607E-01
        1   -1    0    0  503    0  0.00000000000E+00  0.00000000000E+00  1.58591490197E+03  1.58591490197E+03  0.00000000000E+00 0.00000000000E+00  1.
       -1   -1    0    0    0  503  0.00000000000E+00  0.00000000000E+00 -8.99084923758E+00  8.99084923758E+00  0.00000000000E+00 0.00000000000E+00  1.
       23    2    1    2    0    0  4.31808951699E+01  1.18843550193E+01  8.22005355890E+02  8.28398612649E+02  9.24425698805E+01 0.00000000000E+00  1.
       25    2    1    2    0    0 -4.31808951699E+01 -1.18843550193E+01  7.54918696840E+02  7.66507138556E+02  1.25000508063E+02 0.00000000000E+00  1.
       11    1    3    3    0    0 -1.35803884002E+01 -5.28931958672E+00  5.41360784563E+02  5.41556924907E+02  5.11072900539E-04 0.00000000000E+00  1.
      -11    1    3    3    0    0  5.67612835701E+01  1.71736746060E+01  2.80644571326E+02  2.86841687743E+02  5.11012071458E-04 0.00000000000E+00  1.
       23    2    4    4    0    0 -2.43038338852E+01  5.06442605250E+00  2.48359236741E+02  2.53284239962E+02  4.30612469142E+01 0.00000000000E+00  1.
       23    2    4    4    0    0 -1.88770612847E+01 -1.69487810718E+01  5.06559460099E+02  5.13222898594E+02  7.84324703350E+01 0.00000000000E+00  1.
      -13    1    7    7    0    0 -3.25370809281E+01 -6.79837669312E+00  2.02354268485E+02  2.05066186143E+02  1.05659999991E-01 0.00000000000E+00  1.
       13    1    7    7    0    0  8.23324704291E+00  1.18628027456E+01  4.60049682560E+01  4.82180538193E+01  1.05659999989E-01 0.00000000000E+00  1.
      -13    1    8    8    0    0  4.59433181687E+00 -3.18015647781E+01  4.39027117172E+02  4.40201395027E+02  1.05659999655E-01 0.00000000000E+00  1.
       13    1    8    8    0    0 -2.34713931016E+01  1.48527837063E+01  6.75323429266E+01  7.30215035668E+01  1.05660000010E-01 0.00000000000E+00  1.
</event>
  """
  event3 = """
<event>
11  60   1.0000000E+00   1.2500000E+02   7.8125000E-03   1.2380607E-01
        1   -1    0    0  501    0  0.00000000000E+00  0.00000000000E+00  8.38349783822E+01  8.38349783822E+01  0.00000000000E+00 0.00000000000E+00  1.
        2   -1    0    0  502    0  0.00000000000E+00  0.00000000000E+00 -8.69647303563E+02  8.69647303563E+02  0.00000000000E+00 0.00000000000E+00  1.
        4    1    1    2  501    0  4.93534233194E+01 -7.45486758049E+00  2.54822242213E+01  5.60417629563E+01  0.00000000000E+00 0.00000000000E+00  1.
        1    1    1    2  502    0 -4.29482465415E+01  4.39907893858E+01 -7.51475061906E+02  7.53985749267E+02  0.00000000000E+00 0.00000000000E+00  1.
       25    2    1    2    0    0 -6.40517677787E+00 -3.65359218053E+01 -5.98194874970E+01  1.43454769722E+02  1.25000000000E+02 0.00000000000E+00  1.
       23    2    5    5    0    0 -1.61638014503E+01 -3.55963825472E+01 -2.51394501445E+01  1.03431837860E+02  9.24001201399E+01 0.00000000000E+00  1.
       23    2    5    5    0    0  9.75862467247E+00 -9.39539258134E-01 -3.46800373525E+01  4.00229318615E+01  1.74073718437E+01 0.00000000000E+00  1.
      -11    1    6    6    0    0  3.37109433312E+01 -2.97615359833E+01  4.38251799494E+00  4.51816687231E+01  5.11000134768E-04 0.00000000000E+00  1.
       11    1    6    6    0    0 -4.98747447816E+01 -5.83484656388E+00 -2.95219681394E+01  5.82501691374E+01  5.11001208360E-04 0.00000000000E+00  1.
      -13    1    7    7    0    0  1.46596263059E+01  5.33582780943E-01 -2.31337995488E+01  2.73929406894E+01  1.05660000000E-01 0.00000000000E+00  1.
       13    1    7    7    0    0 -4.90100163341E+00 -1.47312203908E+00 -1.15462378037E+01  1.26299911721E+01  1.05660000000E-01 0.00000000000E+00  1.
</event>
  """

  daughters = """
  11 -71.89077865749999319 30.50307494750000004 -47.20025487019999844 91.25012710839999386
  -11 -25.13451734110000046 -18.85931656560000036 -81.42283896300000379 87.27597887359999618
  11 -51.80274100940000181 1.64269040236999997 -41.79162596869999646 66.57899375339999892
  -11 -93.72924763700000028 39.45060783929999815 -92.98363978320000456 137.79506373300000632
  """
  associated = """
  -11 211.33318543799998679 -14.90577872979999974 3.74371777679000006 211.89127619999999297
  12 31.22409920730000010 -37.83127789369999761 1.23465418111000003 49.06805813689999951
  """
  mothers = """
  -1 0.00000000000000000 0.00000000000000000 192.71975508899998886 192.71975508899998886
  2 0.00000000000000000 0.00000000000000000 -451.13974271600000066 451.13974271600000066
  """

  m.setInputEvent(
                  SimpleParticleCollection_t(line.split() for line in daughters.split("\n") if line.split()),
                  SimpleParticleCollection_t(line.split() for line in associated.split("\n") if line.split()),
                  SimpleParticleCollection_t(line.split() for line in mothers.split("\n") if line.split()),
                  True,
                 )
  #or:
  #m.setInputEvent_fromLHE(event1)

  couplings = (
               (1, 0, 0, 0),
               (0, 1, 0, 0),
               (0, 0, 1, 0),
               (0, 0, 0, 1),
               (1, 1.663195, 0, 0),
               (1, 0, 2.55502, 0),
               (1, 0, 0, -12110.20),
              )
  for _ in couplings:
    m.ghz1, m.ghz2, m.ghz4, m.ghz1_prime2 = _
    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Lep_WH)
    prod = m.computeProdP(False)
    m.ghz1, m.ghz2, m.ghz4, m.ghz1_prime2 = _
    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
    dec = m.computeP(False)
    print prod, dec, prod*dec

  print m.computeDecayAngles()
  print "propagator:"
  print "   BW:", m.getXPropagator(TVar.FixedWidth)
  print "  CPS:", m.getXPropagator(TVar.CPS)
