import abc
import re

import ROOT

from mela import Mela, SimpleParticleCollection_t

if __name__ == "__main__":
  import argparse, itertools, sys, unittest
  from mela import TVar
  parser = argparse.ArgumentParser()
  parser.add_argument('--lhefile-hwithdecay')
  parser.add_argument('--lhefile-jhugenvbfvh')
  parser.add_argument('--lhefile-jhugentth')
  parser.add_argument('unittest_args', nargs='*')
  args = parser.parse_args()

class LHEFileBase(object):
  """
  Simple class to iterate through an LHE file and calculate probabilities for each event
  Example usage:
  h1 = ROOT.TH1F("costheta1", "costheta1", 100, -1, 1)
  h2 = ROOT.TH1F("D_0minus", "D_0minus", 100, 0, 1)
  with LHEFile("filename.lhe") as f:
    for event in f:  #event becomes the mela object
      h1.Fill(event.computeDecayAngles().costheta1)
      event.ghz1 = 1
      p0plus = event.computeP()
      event.ghz4 = 1
      p0minus = event.computeP()
      h2.Fill(p0plus / (p0plus + p0minus))
  """
  __metaclass__ = abc.ABCMeta

  __melas = {}

  def __init__(self, filename, *melaargs, **kwargs):
    self.isgen = kwargs.pop("isgen", True)
    reusemela = kwargs.pop("reusemela", False)
    if kwargs: raise ValueError("Unknown kwargs: " + ", ".join(kwargs))
    self.filename = filename
    if reusemela and melaargs in self.__melas:
      self.mela = self.__melas[melaargs]
    else:
      self.__melas[melaargs] = self.mela = Mela(*melaargs)
    self.f = open(self.filename)
  def __enter__(self, *args, **kwargs):
    self.f.__enter__(*args, **kwargs)
    return self
  def __exit__(self, *args, **kwargs):
    return self.f.__exit__(*args, **kwargs)

  def __iter__(self):
    event = ""
    for linenumber, line in enumerate(self.f, start=1):
      if "<event>" not in line and not event:
        continue
      event += line
      if "</event>" in line:
        try:
          self._setInputEvent(event)
          yield self
          event = ""
        except GeneratorExit:
          raise
        except:
          print "On line", linenumber
          raise
        finally:
          try:
            self.mela.resetInputEvent()
          except:
            pass

  @abc.abstractmethod
  def _setInputEvent(self, event): pass

  @classmethod
  def _LHEclassattributes(cls):
    return "filename", "f", "mela", "isgen"

  def __getattr__(self, attr):
    if attr == "mela": raise RuntimeError("Something is wrong, trying to access mela before it's created")
    return getattr(self.mela, attr)
  def __setattr__(self, attr, value):
    if attr in self._LHEclassattributes():
      super(LHEFileBase, self).__setattr__(attr, value)
    else:
      setattr(self.mela, attr, value)

class LHEFile_Hwithdecay(LHEFileBase):
  def _setInputEvent(self, event):
    self.mela.setInputEvent_fromLHE(event, self.isgen)

class LHEFile_JHUGenVBFVH(LHEFileBase):
  def _setInputEvent(self, event):
    self.ids = [int(line.split()[0]) for line in event.strip().split("\n")[2:-1]]
    if not self.isgen:
      event = re.sub("^( *-?)([12345]|21) ", r"\g<1>0 ", event, flags=re.MULTILINE)
    particles = event.strip().split("\n")[2:-1]  #remove the first 2 lines, <event> and event info, and the last line, </event>
    self.daughters = SimpleParticleCollection_t([particle for particle in particles if int(particle.split()[0]) == 25])
    self.associated = SimpleParticleCollection_t([particle for particle in particles
                             if abs(int(particle.split()[0])) in (0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16)
                             and int(particle.split()[1]) == 1])
    self.mothers = SimpleParticleCollection_t([particle for particle in particles if int(particle.split()[1]) == -1])
    assert len(self.daughters) == 1 and len(self.associated) == len(self.mothers) == 2
    self.mela.setInputEvent(self.daughters, self.associated, self.mothers, self.isgen)
  @classmethod
  def _LHEclassattributes(cls):
    return super(LHEFile_JHUGenVBFVH, cls)._LHEclassattributes() + ("daughters", "mothers", "associated", "ids")

class LHEFile_JHUGenttH(LHEFileBase):
  def _setInputEvent(self, event):
    self.ids = [int(line.split()[0]) for line in event.strip().split("\n")[2:-1]]
    if not self.isgen:
      event = re.sub("^( *-?)([12345]|21) ", r"\g<1>0 ", event, flags=re.MULTILINE)
    particles = event.strip().split("\n")[2:-1]  #remove the first 2 lines, <event> and event info, and the last line, </event>
    self.daughters = SimpleParticleCollection_t([particle for particle in particles if int(particle.split()[0]) == 25])
    self.associated = SimpleParticleCollection_t([particle for particle in particles
                             if abs(int(particle.split()[0])) in (0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16)
                             and int(particle.split()[1]) == 1])
    self.mothers = SimpleParticleCollection_t([particle for particle in particles if int(particle.split()[1]) == -1])
    assert len(self.daughters) == 1 and len(self.associated) == 6 and len(self.mothers) == 2
    self.mela.setInputEvent(self.daughters, self.associated, self.mothers, self.isgen)
  @classmethod
  def _LHEclassattributes(cls):
    return super(LHEFile_JHUGenttH, cls)._LHEclassattributes() + ("daughters", "mothers", "associated", "ids")


if __name__ == '__main__':
  class TestLHEFiles(unittest.TestCase):
    @unittest.skipUnless(args.lhefile_hwithdecay, "needs --lhefile-hwithdecay argument")
    def testHwithDecay(self):
      with LHEFile_Hwithdecay(args.lhefile_hwithdecay) as f:
        for event, i in itertools.izip(f, range(10)):
          event.ghz1 = 1
          event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
          prob = event.computeP()
          self.assertNotEqual(prob, 0)
          print prob, event.computeDecayAngles()

    @unittest.skipUnless(args.lhefile_jhugenvbfvh, "needs --lhefile-jhugenvbfvh argument")
    def testJHUGenVBFVH(self):
      with LHEFile_JHUGenVBFVH(args.lhefile_jhugenvbfvh) as f:
        for event, i in itertools.izip(f, range(10)):
          event.ghz1 = 1
          if any(11 <= abs(p.first) <= 16 for p in event.associated):
            if sum(p.first for p in event.associated) == 0:
              event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Lep_ZH)
            else:
              event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.Lep_WH)
          else:
            event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
          prob = event.computeP()
          self.assertNotEqual(prob, 0)
          print prob, event.computeVBFAngles(), event.computeVHAngles()

    @unittest.skipUnless(args.lhefile_jhugentth, "needs --lhefile-jhugentth argument")
    def testJHUGenttH(self):
      with LHEFile_JHUGenttH(args.lhefile_jhugentth) as f:
        for event, i in itertools.izip(f, range(10)):
          event.ghz1 = 1
          event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
          prob = event.computeP()
          self.assertNotEqual(prob, 0)
          print prob, event.computeDecayAngles()

  unittest.main(argv=[sys.argv[0]]+args.unittest_args)
