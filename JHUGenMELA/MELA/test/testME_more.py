#!/usr/bin/env python

"""
test python and some MEs
"""

import os
import random
import sys
import unittest

from mela import Mela, TVar

event1_ggH = """
<event>
 9   0   1.0000000E+00   6.2499553E+01   7.8125000E-03   1.3828638E-01
       21   -1    0    0  501  502  0.00000000000E+00  0.00000000000E+00  5.37086667033E+02  5.37086667033E+02  0.00000000000E+00 0.00000000000E+00  1.
       21   -1    0    0  502  501  0.00000000000E+00  0.00000000000E+00 -7.27293062837E+00  7.27293062837E+00  0.00000000000E+00 0.00000000000E+00  1.
       25    2    1    2    0    0  2.66453525910E-15  0.00000000000E+00  5.29813736405E+02  5.44359597662E+02  1.24999105129E+02 0.00000000000E+00  1.
       23    2    3    3    0    0 -1.07287063341E+00  2.03413884495E+01  2.12763358296E+02  2.18095469074E+02  4.33873698414E+01 0.00000000000E+00  1.
       23    2    3    3    0    0  1.07287063341E+00 -2.03413884495E+01  3.17050378109E+02  3.26264128588E+02  7.42456477420E+01 0.00000000000E+00  1.
       11    1    4    4    0    0 -4.96385765951E+00  1.62933935971E+01  2.11696003997E+02  2.12380113632E+02  5.10998597706E-04 0.00000000000E+00  1.
      -11    1    4    4    0    0  3.89098702611E+00  4.04799485238E+00  1.06735429842E+00  5.71535544142E+00  5.10999912595E-04 0.00000000000E+00  1.
       12    1    5    5    0    0  8.04215953543E+00 -2.29213235677E+01  2.77968427622E+01  3.69151442044E+01  4.63538696652E-06 0.00000000000E+00  1.
      -12    1    5    5    0    0 -6.96928890202E+00  2.57993511828E+00  2.89253535347E+02  2.89348984384E+02  5.39479660939E-06 0.00000000000E+00  1.
</event>
"""

event2_VBF = """
<event>
11  60   1.0000000E+00   1.2500000E+02   7.8125000E-03   1.2380607E-01
       -4   -1    0    0  501    0  0.00000000000E+00  0.00000000000E+00  3.03951895455E+01  3.03951895455E+01 0.00000000000E+00 0.00000000000E+00  1.
        1   -1    0    0  502    0  0.00000000000E+00  0.00000000000E+00 -7.16454865818E+02  7.16454865818E+02 0.00000000000E+00 0.00000000000E+00  1.
       -4    1    1    2  501    0 -6.29623305954E+01 -4.67693451618E+01 -5.85204921720E+01  9.78584422772E+01 0.00000000000E+00 0.00000000000E+00  1.
        1    1    1    2  502    0  9.29867095677E+00  2.61406963578E+01 -9.53450713739E+01  9.92999694004E+01 0.00000000000E+00 0.00000000000E+00  1.
       25    2    1    2    0    0  5.36636596386E+01  2.06286488040E+01 -5.32194112726E+02  5.49691643686E+02 1.25000000000E+02 0.00000000000E+00  1.
       23    2    5    5    0    0 -4.38443764422E-01 -1.86323403603E+01 -9.56577645198E+01  9.86840261342E+01 1.55168540914E+01 0.00000000000E+00  1.
       23    2    5    5    0    0  5.41021034030E+01  3.92609891643E+01 -4.36536348206E+02  4.51007617552E+02 9.15173476546E+01 0.00000000000E+00  1.
      -11    1    6    6    0    0  3.31752759767E+00  1.23639112123E+00 -8.43121325514E+00  9.14439771558E+00 5.10999986247E-04 0.00000000000E+00  1.
       11    1    6    6    0    0 -3.75597136210E+00 -1.98687314816E+01 -8.72265512647E+01  8.95396284186E+01 5.11001208360E-04 0.00000000000E+00  1.
      -13    1    7    7    0    0  6.93884009242E+01  3.00845351746E+01 -2.00352143436E+02  2.14151399675E+02 1.05660000076E-01 0.00000000000E+00  1.
       13    1    7    7    0    0 -1.52862975212E+01  9.17645398978E+00 -2.36184204771E+02  2.36856217877E+02 1.05659999949E-01 0.00000000000E+00  1.
</event>
"""

event3_ZH = """
<event>
12  50   1.0000000E+00   1.2500000E+02   7.8125000E-03   1.2380607E-01
        1   -1    0    0  503    0  0.00000000000E+00  0.00000000000E+00  5.98124082110E+01  5.98124082110E+01  0.00000000000E+00 0.00000000000E+00  1.
       -1   -1    0    0    0  503  0.00000000000E+00  0.00000000000E+00 -3.44771438947E+02  3.44771438947E+02  0.00000000000E+00 0.00000000000E+00  1.
       23    2    1    2    0    0 -1.33768582227E+01  8.88541882213E+01 -1.21353594052E+02  1.79939405232E+02  9.78646395061E+01 0.00000000000E+00  1.
       25    2    1    2    0    0  1.33768582227E+01 -8.88541882213E+01 -1.63605436683E+02  2.24644441925E+02  1.24997517076E+02 0.00000000000E+00  1.
       13    1    3    3    0    0  3.17738029751E+01 -1.54337455361E+00 -2.87056115717E+01  4.28483355599E+01  1.05660000003E-01 0.00000000000E+00  1.
      -13    1    3    3    0    0 -4.51506611978E+01  9.03975627749E+01 -9.26479824805E+01  1.37091069673E+02  1.05660000012E-01 0.00000000000E+00  1.
       23    2    4    4    0    0  4.40015216426E+01 -3.16198925284E+01 -5.72188020543E+01  7.95131988002E+01  1.06021679173E+01 0.00000000000E+00  1.
       23    2    4    4    0    0 -3.06246634199E+01 -5.72342956929E+01 -1.06386634629E+02  1.45131243125E+02  7.43728921737E+01 0.00000000000E+00  1.
      -11    1    8    8    0    0  1.78742924427E+01 -1.12102312034E+01 -2.03424438331E+00  2.11966451221E+01  5.10999863718E-04 0.00000000000E+00  1.
       11    1    7    7    0    0  3.86248305415E+01 -3.12040320381E+01 -4.96421050079E+01  7.02133017289E+01  5.11000665199E-04 0.00000000000E+00  1.
      -11    1    7    7    0    0  5.37669110106E+00 -4.15860490296E-01 -7.57669704638E+00  9.29989707133E+00  5.11000011708E-04 0.00000000000E+00  1.
       11    1    8    8    0    0 -4.84989558626E+01 -4.60240644895E+01 -1.04352390245E+02  1.23934598003E+02  5.10996863057E-04 0.00000000000E+00  1.
</event>
"""

class TestMela(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.m = Mela()
  def setUp(self):
    """the seed will be different every time, but reproducible"""
    try:
        random.seed(os.environ["SEED"] + os.environ["CIRCLE_BUILD_NUM"])
    except KeyError:
        random.seed(34567)
  def tearDown(self):
    self.m.resetInputEvent()
  def testcontact_decay(self):
    self.runcontact(event1_ggH, False, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
  def testcontact_VBF(self):
    self.runcontact(event2_VBF, True, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
  def testcontact_ZH(self):
    self.runcontact(event3_ZH, True, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Lep_ZH)
  def runcontact(self, event, isprod, *setprocessargs):
    """
    correspondence between a1 L1 L1Zg and contact terms
    """

    M_Z = 91.1876
    Ga_Z = 2.4952
    aL_lep = -0.53762
    aR_lep = 0.46238
    aL_nu = 1
    aR_nu = 0

    sitW2 = 0.23119
    aL_up = (-2 * sitW2 * 2./3 + 1)
    aR_up = -2 * sitW2 * 2./3
    aL_dn = (-2 * sitW2 * -1./3 - 1)
    aR_dn = -2 * sitW2 * -1./3

    e = 0.8431872482432357  # = cL_lep = cR_lep from mod_Parameters
    c_up = -2./3 * e
    c_dn = 1./3 * e
    L1 = 10000.

    m = self.m
    computeP = m.computeProdP if isprod else m.computeP

    ghz1, ghz1_prime2, ghzgs1_prime2 = random.uniform(-1, 1), random.uniform(-10000, 10000), random.uniform(-10000, 10000)
    m.setInputEvent_fromLHE(event, True)

    m.setProcess(*setprocessargs)
    m.ghz1 = ghz1
    m.ghz1_prime2 = ghz1_prime2
    m.ghzgs1_prime2 = ghzgs1_prime2
    me1 = computeP()

    m.setProcess(*setprocessargs)
    m.ghz1 = ghz1 + 2 * ghz1_prime2 * (M_Z**2 - 1j*M_Z*Ga_Z)/L1**2
    m.ghzzp1 = M_Z**2/L1**2
    m.ezp_El_left = m.ezp_Mu_left = m.ezp_Ta_left = aL_lep * ghz1_prime2 + e * ghzgs1_prime2
    m.ezp_El_right = m.ezp_Mu_right = m.ezp_Ta_right = aR_lep * ghz1_prime2 + e * ghzgs1_prime2
    m.ezp_Up_left = m.ezp_Chm_left = m.ezp_Top_left = aL_up * ghz1_prime2 + c_up * ghzgs1_prime2
    m.ezp_Up_right = m.ezp_Chm_right = m.ezp_Top_right = aR_up * ghz1_prime2 + c_up * ghzgs1_prime2
    m.ezp_Dn_left = m.ezp_Str_left = m.ezp_Bot_left = aL_dn * ghz1_prime2 + c_dn * ghzgs1_prime2
    m.ezp_Dn_right = m.ezp_Str_right = m.ezp_Bot_right = aR_dn * ghz1_prime2 + c_dn * ghzgs1_prime2
    m.ezp_NuE_left = aL_nu * ghz1_prime2
    m.ezp_NuE_right = aR_nu * ghz1_prime2
    me2 = computeP()

    self.assertEquals(me1, me2)
    self.assertNotEquals(me1, 0)

  def testzp_decay(self):
    self.runVprime(event1_ggH, False, False, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
  def testzp_VBF(self):
    self.runVprime(event2_VBF, True, False, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
  def testzp_ZH(self):
    self.runVprime(event3_ZH, True, False, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Lep_ZH)
  def testzpzp_decay(self):
    self.runVprime(event1_ggH, False, True, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
  def testzpzp_VBF(self):
    self.runVprime(event2_VBF, True, True, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
  def testzpzp_ZH(self):
    self.runVprime(event3_ZH, True, True, TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Lep_ZH)
  def runVprime(self, event, isprod, usevpvp, *setprocessargs):
    """
    correspondence between a1 L1 L1Zg and contact terms
    """

    M_Z = 91.1876
    Ga_Z = 2.4952
    M_W = 80.399
    Ga_W = 2.085
    sitW2 = 0.23119

    aL_lep = -0.53762
    aR_lep = 0.46238
    aL_nu = 1
    aR_nu = 0

    sitW2 = 0.23119
    aL_up = (-2 * sitW2 * 2./3 + 1)
    aR_up = -2 * sitW2 * 2./3
    aL_dn = (-2 * sitW2 * -1./3 - 1)
    aR_dn = -2 * sitW2 * -1./3

    bL = (2*(1-sitW2))**.5

    m = self.m
    computeP = m.computeProdP if isprod else m.computeP
    m.setInputEvent_fromLHE(event, True)

    m.setProcess(*setprocessargs)
    m.ghz1 = 1
    me1 = computeP()

    m.setProcess(*setprocessargs)

    if usevpvp:
        m.ghzpzp1 = 1
    else:
        m.ghzzp1 = 0.5
    m.M_Zprime = M_Z
    m.Ga_Zprime = Ga_Z
    m.M_Wprime = M_W
    m.Ga_Wprime = Ga_W

    m.ezp_El_left = m.ezp_Mu_left = m.ezp_Ta_left = aL_lep
    m.ezp_El_right = m.ezp_Mu_right = m.ezp_Ta_right = aR_lep
    m.ezp_Up_left = m.ezp_Chm_left = m.ezp_Top_left = aL_up
    m.ezp_Up_right = m.ezp_Chm_right = m.ezp_Top_right = aR_up
    m.ezp_Dn_left = m.ezp_Str_left = m.ezp_Bot_left = aL_dn
    m.ezp_Dn_right = m.ezp_Str_right = m.ezp_Bot_right = aR_dn
    m.ezp_NuE_left = aL_nu
    m.ezp_NuE_right = aR_nu

    m.ewp_El_left = m.ewp_Mu_left = m.ewp_Ta_left = m.ewp_Up_left =  m.ewp_Chm_left = m.ewp_Top_left = bL

    me2 = computeP()

    self.assertEquals(me1, me2)
    self.assertNotEquals(me1, 0)

if __name__ == "__main__":
  unittest.main()
