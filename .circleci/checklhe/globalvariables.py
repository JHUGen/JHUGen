import particlecategory
import particle
import config
import collections
import cPickle

eventcounter = collections.Counter()

startedinit = False
finishedinit = False

class GlobalVariablesMinimal:
    #for stuff that's needed before GlobalVariables is inited
    #I'm aware this is ugly
    def __init__(self):
        class ParticleCategory:
            def __init__(self, particles):
                self.particles = particles
            def ids(self):
                return self.particles

        self.neutralbosons = ParticleCategory([21, 22, 23, 25, 32, 39])

class GlobalVariables:
    def __init__(self):
        self.startedinit = False
        self.finishedinit = False

    #Not __init__, I want the object to exist even before the variables are assigned
    def init(self):
        if self.startedinit:
            return
        self.startedinit = True

        #Particle categories
        self.neutralbosons = particlecategory.ParticleCategory([21, 22, 23, 25, 32, 39, 625], Csymmetric = False)

        self.electrons = particlecategory.ParticleCategory([11])
        self.muons = particlecategory.ParticleCategory([13])
        self.taus = particlecategory.ParticleCategory([15])
        self.emu = self.electrons.union(self.muons)
        self.leptons = self.emu.union(self.taus)

        self.neutrinos = particlecategory.ParticleCategory([12, 14, 16])

        self.down    = particlecategory.ParticleCategory([1])
        self.up      = particlecategory.ParticleCategory([2])
        self.strange = particlecategory.ParticleCategory([3])
        self.charm   = particlecategory.ParticleCategory([4])
        self.bottom  = particlecategory.ParticleCategory([5])
        self.top     = particlecategory.ParticleCategory([6])
        self.uptypequarks = self.up.union(self.charm).union(self.top)
        self.downtypequarks = self.down.union(self.strange).union(self.bottom)
        self.quarks = self.uptypequarks.union(self.downtypequarks)

        self.gluon = particlecategory.ParticleCategory([21])
        self.photon = particlecategory.ParticleCategory([22])
        self.weakbosons = particlecategory.ParticleCategory([23, 24])
        self.Z = particlecategory.ParticleCategory([23])
        self.W = particlecategory.ParticleCategory([24])
        self.higgs = particlecategory.ParticleCategory([25, 32, 39])

        self.jets = self.quarks.difference(self.top).union(self.gluon)
        print "initialized particle categories"

        config.init()
        self.finishedinit = True

def init():
    global globalvariables, startedinit, finishedinit
    if startedinit:
        return
    startedinit = True

    globalvariables = GlobalVariables()
    globalvariables.init()

    config.init()
    finishedinit = True

