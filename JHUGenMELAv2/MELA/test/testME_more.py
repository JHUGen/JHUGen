#!/usr/bin/env python

"""
test python and some MEs
"""

import os
import random
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
       13    1    5    5    0    0  8.04215953543E+00 -2.29213235677E+01  2.77968427622E+01  3.69151442044E+01  4.63538696652E-06 0.00000000000E+00  1.
      -13    1    5    5    0    0 -6.96928890202E+00  2.57993511828E+00  2.89253535347E+02  2.89348984384E+02  5.39479660939E-06 0.00000000000E+00  1.
</event>
"""

class TestMela(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.m = Mela()
  def setUp(self):
    """the seed will be different every time, but reproducible"""
    random.seed(os.environ["SEED"] + os.environ["CIRCLE_BUILD_NUM"])
  def tearDown(self):
    self.m.resetInputEvent()
  def testcontact(self):
    """
    correspondence between a1 L1 L1Zg and contact terms
    """

    M_Z = 91.1876
    Ga_Z = 2.4952
    aL_lep = -0.53762
    aR_lep = 0.46238
    aL_up = -2 * 0.23119 * 2./3
    aR_up = -2 * 0.23119 * 2./3 + 1
    aL_dn = -2 * 0.23119 * -1./3
    aR_dn = -2 * 0.23119 * -1./3 - 1
    e = 0.8431872482432357  # = cL_lep = cR_lep from mod_Parameters
    L1 = 10000.

    m = self.m

    ghz1, ghz1_prime2, ghzgs1_prime2 = random.uniform(-1, 1), random.uniform(-10000, 10000), random.uniform(-10000, 10000)
    m.setInputEvent_fromLHE(event1_ggH, True)

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
    m.ghz1 = ghz1
    m.ghz1_prime2 = ghz1_prime2
    m.ghzgs1_prime2 = ghzgs1_prime2
    me1 = m.computeP()

    m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
    m.ghz1 = ghz1 + 2 * ghz1_prime2 * (M_Z**2 - 1j*M_Z*Ga_Z)/L1**2
    m.ghzzp1 = M_Z**2/L1**2
    m.ezp_El_left = m.ezp_Mu_left = m.ezp_Ta_left = aL_lep * ghz1_prime2 + e * ghzgs1_prime2
    m.ezp_El_right = m.ezp_Mu_right = m.ezp_Ta_right = aR_lep * ghz1_prime2 + e * ghzgs1_prime2
    m.ezp_Up_left = m.ezp_Chm_left = m.ezp_Top_left = aL_up * ghz1_prime2 + e*2/3 * ghzgs1_prime2
    m.ezp_Up_right = m.ezp_Chm_right = m.ezp_Top_right = aR_up * ghz1_prime2 + e*2/3 * ghzgs1_prime2
    m.ezp_Dn_left = m.ezp_Str_left = m.ezp_Bot_left = aL_dn * ghz1_prime2 + e*-1/3 * ghzgs1_prime2
    m.ezp_Dn_right = m.ezp_Str_right = m.ezp_Bot_right = aR_dn * ghz1_prime2 + e*-1/3 * ghzgs1_prime2
    me2 = m.computeP()

    self.assertEquals(me1, me2)

if __name__ == "__main__":
  unittest.main()
