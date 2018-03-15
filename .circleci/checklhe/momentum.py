import usefulstuff
import config
import ROOT
from math import acos, atan2, copysign, log, tan

class Momentum(ROOT.TLorentzVector):
    def __init__(self, event, *args, **kwargs):
        ROOT.TLorentzVector.__init__(self, *args, **kwargs)
        for component in self.Px(), self.Py(), self.Pz(), self.E():
            if not usefulstuff.isfinite(component):
                raise ValueError
        self.event = event
        if event is not None:
            event.momenta.append(self)
    def momentum(self):
        return self        #so that particle and momentum can be used symmetrically
    def invmass(self):
        return self.M()
    def __neg__(self):
        return Momentum(self.event, -self.Px(), -self.Py(), -self.Pz(), -self.E())
    def __add__(self, other):
        return Momentum(self.event, self.Px()+other.Px(), self.Py()+other.Py(), self.Pz()+other.Pz(), self.E()+other.E())
    def __sub__(self, other):
        return self + (-other)
    def __eq__(self, other):
        return (self-other).P()**2 + (self-other).E()**2 < config.momentumtolerance**2
    def __ne__(self, other):
        return not (self == other)
    def __str__(self):
        return "(%s, %s, %s, %s)" % (self.Px(), self.Py(), self.Pz(), self.E())

class Frame(object):
    def __init__(self, event):
        self.x = Momentum(event, 1, 0, 0, 0)
        self.y = Momentum(event, 0, 1, 0, 0)
        self.z = Momentum(event, 0, 0, 1, 0)
        self.t = Momentum(event, 0, 0, 0, 1)

    def Boost(self, *args):
        for v in (self.x, self.y, self.z, self.t):
            v.Boost(*args)
    def Rotate(self, *args):
        for v in (self.x, self.y, self.z, self.t):
            v.Rotate(*args)
