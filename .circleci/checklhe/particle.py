from fractions import Fraction
from collections import Counter
import momentum
import usefulstuff
import itertools
import particletype
import particlecategory
import globalvariables
import config
import vertex
import color

class BadParticleLineError(Exception):
    pass

class Particle(particletype.ParticleType):
    def __init__(self, line, ev):
        self.__line = line
        self.ev = ev
        self.miscellaneouschecks = []
        data = line.split()
        if len(data) < 13:
            raise BadParticleLineError("This line is not a valid particle line.  Not enough entries.")

        try:
            super(Particle, self).__init__(int(data[0]))
        except ValueError:
            raise BadParticleLineError
        try:
            self.__status = int(data[1])
        except ValueError:
            self.__status = -999
            self.miscellaneouschecks.append("Status for " + str(self) + " is " + data[1] + " instead of a number!")

        try:
            self.__color = int(data[4])
        except ValueError:
            self.miscellaneouschecks.append("Color for " + str(self) + " is " + data[4] + " instead of a number!")
            self.__color = 0
        try:
            self.__anticolor = int(data[5])
        except ValueError:
            self.miscellaneouschecks.append("Anticolor for " + str(self) + " is " + data[5] + " instead of a number!")
            self.__color = 0

        try:
            self.__momentum = momentum.Momentum(ev, float(data[6]), float(data[7]), float(data[8]), float(data[9]))
        except ValueError:
            self.miscellaneouschecks.append("Momentum for " + str(self) + " is " + "(" + ", ".join(data[6:10]) + ")!")
            self.__momentum = momentum.Momentum(ev, 0, 0, 0, 0)

        try:
            self.__lhemass = float(data[10])
            if not usefulstuff.isfinite(self.__lhemass): raise ValueError
        except ValueError:
            self.miscellaneouschecks.append("Mass for " + str(self) + " is " + data[11] + " instead of a number!")
            self.__lhemass = 0
        try:
            self.__lifetime = float(data[11])
            if not usefulstuff.isfinite(self.__lifetime): raise ValueError
        except ValueError:
            self.miscellaneouschecks.append("Lifetime for " + str(self) + " is " + data[11] + " instead of a number!")
            self.__lifetime = 0
        try:
            self.__spin = float(data[12])
            if not usefulstuff.isfinite(self.__spin): raise ValueError
        except ValueError:
            self.miscellaneouschecks.append("Spin for " + str(self) + " is " + data[12] + " instead of a number!")
            self.__spin = 9

        self.__mothers = usefulstuff.printablefrozenset([int(data[2]), int(data[3])])
        self.__mothersset = False
        self.__kids = usefulstuff.printablelist([])

        ev.particlelist.append(self)

    def setmothers(self):
        if self.__mothersset:
            return
        self.__mothersset = True
        self.__mothers = usefulstuff.printablefrozenset([self.ev.particlelist[m] for m in self.__mothers])
        self.__mothers.setmothers()
        for particle in self.ev.particlelist[1:]:
            if particle in self.mothers():
                particle.kids().append(self)
        self.startvertex = self.ev.vertices[self.mothers()]
        self.endvertex = None
        for mother in self.mothers():
            if mother is None:
                continue
            if mother.endvertex is None:
                mother.endvertex = self.startvertex
            elif mother.endvertex is not self.startvertex:
                self.miscellaneouschecks.append(str(mother) + " decays in multiple vertices!")
        if self.startvertex is not None:
            self.startvertex.addkid(self)

        if self.color() != 0:
            self.ev.colors[self.color()].addparticle(self)
        if self.anticolor() != 0:
            self.ev.colors[self.anticolor()].addparticle(self)

    def __eq__(self, other):
       return self is other
    def __ne__(self, other):
       return not self == other
    def __hash__(self):
       return hash(self.__line)

    def status(self):
        return self.__status
    def color(self):
        return self.__color
    def anticolor(self):
        return self.__anticolor
    def mothers(self):
        return self.__mothers
    def kids(self):
        return self.__kids
    def momentum(self):
        return self.__momentum
    def invmass(self):
        return self.__momentum.M()
    def lhemass(self):
        return self.__lhemass

    def Px(self):
        return self.__momentum.Px()
    def Py(self):
        return self.__momentum.Py()
    def Pz(self):
        return self.__momentum.Pz()
    def E(self):
        return self.__momentum.E()
    def Pt(self):
        return self.__momentum.Pt()
    def Eta(self):
        return self.__momentum.Eta()
    def CosTheta(self):
        return self.__momentum.CosTheta()
    def Phi(self):
        return self.__momentum.Phi()
    def Rapidity(self):
        return self.__momentum.Rapidity()

    def Boost(self, *args):
        self.__momentum.Boost(*args)
    def Rotate(self, *args):
        self.__momentum.Rotate(*args)
    def Vect(self):
        return self.__momentum.Vect()
    def BoostVector(self):
        return self.__momentum.BoostVector()

class ParticleCounter(Counter):
    def __init__(self, particles):
        try:                    #if particles is another ParticleCounter
            super(ParticleCounter, self).__init__([particletype.ParticleType(p) for p in particles.elements() if p is not None])
        except AttributeError:
            try:                #if it's a list
                super(ParticleCounter, self).__init__([particletype.ParticleType(p) for p in particles if p is not None])
            except TypeError:   #if it's a single particle
                super(ParticleCounter, self).__init__([particletype.ParticleType(particles)])
    def count(self, whattocount):
        return len([p for p in self.elements() if p in whattocount])
    def charge(self):
        return sum([p.charge() for p in self.elements()])
    def baryonnumber(self):
        return sum([p.baryonnumber() for p in self.elements()])
    def leptonnumber(self, generation):
        return sum([p.leptonnumber(generation) for p in self.elements()])
    def __hash__(self):
        return hash(tuple(sorted(hash(p) for p in self.elements())))
    def __str__(self):
        return " ".join(str(p) for p in self.elements())
    def __eq__(self, other):
        return hash(self) == hash(other)
    def __ne__(self, other):
        return not self == other

class DecayType(ParticleCounter):
    def __init__(self, particle, level = None):
        try:
            decayparticles = particle.particles   #see if it's itself a decaytype
        except AttributeError:
            try:
                decayparticles = [p for p in particle if p is not None]
            except TypeError:
                decayparticles = [particle]
        done = False
        i = 0
        while not done and (level is None or i < level):
            i += 1
            done = True
            for p in decayparticles[:]:
                if p.kids():
                    done = False
                    decayparticles.remove(p)
                    decayparticles += p.kids()
        self.particles = usefulstuff.printablelist(decayparticles)
        super(DecayType, self).__init__(decayparticles)

class EventCount(object):
    def __init__(self, name = "", subcategories = None, dontprintifnonew = False, dontprintifparentisnt = False):
        self.name = name
        self.__active = True
        self.dontprintifnonew = dontprintifnonew
        self.dontprintifparentisnt = dontprintifparentisnt
        self.subcategories = []
        if subcategories is not None:
            self.subcategories = subcategories

    def __hash__(self):
        return hash(self.name)

    def getcount(self, eventcounter):
        return eventcounter[self]

    def printcount(self, eventcounter):
        count = self.getcount(eventcounter)
        total = globalvariables.globalvariables.anyevent.getcount(eventcounter)
        toreactivate = []
        if self.dontprintifnonew:
            for subcategory in self.subcategories:
                if self.getcount(eventcounter) == subcategory.getcount(eventcounter):
                    toreactivate.append(self)
                    self.deactivate()
                    break
        if self.__active and count > 0:
            joiner = "\n    "
            result = "%s %s events (%s%%)" % (count, self.name, 100.0*count/total)
            result += joiner
        else:
            result = ""
            joiner = "\n"
            for subcategory in self.subcategories:
                if subcategory.dontprintifparentisnt and subcategory.isactive():
                    toreactivate.append(subcategory)
                    subcategory.deactivate()

        for subcategory in self.subcategories:
            result += joiner.join([line for line in subcategory.printcount(eventcounter).split("\n")])

        for category in toreactivate:
            category.activate()

        return result.rstrip(" ")

    def increment(self, eventcounter):
        if self.__active:
            eventcounter[self] += 1

    def activate(self, recursive = False):
        self.__active = True
        if recursive:
            for s in self.subcategories:
                s.activate(recursive)
    def deactivate(self, recursive = False):
        self.__active = False
        if recursive:
            for s in self.subcategories:
                s.deactivate(recursive)
    def isactive(self):
        return self.__active

class DecayFamily(EventCount, set):
    def __init__(self, decaytypes, charge = None, baryonnumber = None, leptonnumber = (None, None, None), name = "", subcategories = None, Csymmetric = True):
        finallist = []
        secondarylist = []

        for d in decaytypes:
            try:
                secondarylist.append(ParticleCounter(d))
            except TypeError:
                for tple in itertools.product(*[particlecategory.ParticleCategory(p, Csymmetric) for p in d]):
                    secondarylist.append(ParticleCounter(tple))

        for d in secondarylist:
            if charge is None or charge == d.charge():
                if baryonnumber is None or d.baryonnumber() == baryonnumber:
                    if all(leptonnumber[i] is None or d.leptonnumber(i+1) == leptonnumber[i] for i in range(3)):
                        finallist.append(d)

        set.__init__(self, finallist)
        EventCount.__init__(self, name, subcategories)

    def increment(self, decay, eventcounter):
        incremented = []
        if len(self) == 0 or decay in self:
            super(DecayFamily, self).increment(eventcounter)
            incremented.append(self)
        for subcategory in self.subcategories:
            subincremented = subcategory.increment(decay, eventcounter)
            incremented += subincremented
            #to check
            if (self not in incremented) and subincremented:
                raise RuntimeError("Subcategory '%s' of category '%s' is filled even though its parent is not!" % (subincremented[0].name, self.name))
        return incremented

    def __contains__(self, other):
        if other is None:
            return False
        if super(DecayFamily, self).__contains__(ParticleCounter(other)):
            return True
        try:
            newother = DecayType(other, 1)
        except AttributeError:
            return False
        if newother != ParticleCounter(other):
            return self.__contains__(newother)
        else:
            return False
    def __hash__(self):
        return hash((self.name, tuple(d for d in self)))
    def __str__(self):
        return ";    ".join(str(decay) for decay in self)
