import particle
import particletype
import particlecategory
import globalvariables
import usefulstuff
import config
import momentum
import vertex
import color
import ROOT
from math import copysign, acos

class Event:

#############################
#        Event setup        #
#############################

    def __init__(self, linenumber, processidlist):
        self.linenumber = linenumber
        self.processidlist = processidlist
        self.firstline = None
        self.particlelist = usefulstuff.printablelist([None])
        self.vertices = vertex.Vertices()
        self.colors = color.Colors()
        self.done = False
        self.momenta = []

    def setfirstline(self, firstline):
        if self.firstline is not None:
            raise ValueError("The first line for this event is already set to:\n" + self.firstline + "\nand cannot be set again to:\n")
        self.firstline = firstline

    def addparticle(self, particleline):
        if self.firstline is None:
            raise ValueError("The first line for this event has not been set yet!")
        if self.done:
            raise ValueError("finished() has already been called for this event, so no more particles can be added!")
        particle.Particle(particleline, self)

    def finished(self):
        if self.firstline is None:
            raise ValueError("The first line for this event has not been set yet!")
        self.particlelist.setmothers()
        self.particlelist.remove(None)
        self.particlecounter = particle.ParticleCounter(self.particlelist)

        self.decaylist = [p for p in self.particlelist if p.kids()]
        self.incoming = [p for p in self.particlelist if p.status() == -1]

        self.miscellaneouschecks = sum(([a + " " + str(self.linenumber) for a in b] for b in self.particlelist[1:].miscellaneouschecks),[])

        self.checkfunctions = []
        if config.checkfirstline:
            self.checkfunctions.append(self.checkfirstline)
        if config.checkstatus:
            self.checkfunctions.append(self.checkstatus)
        if config.checkinvmass:
            self.checkfunctions.append(self.checkinvmass)
        if config.checkPDGmass:
            self.checkfunctions.append(self.checkPDGmass)
        if config.checkmomentum:
            self.checkfunctions.append(self.checkmomentum)
        if config.checkcharge:
            self.checkfunctions.append(self.checkcharge)
        if config.checkcolor:
            self.checkfunctions.append(self.checkcolor)

    def count(self, whattocount):
        return self.particlecounter.count(whattocount)

    def check(self):
        checks = self.miscellaneouschecks + [chk() for chk in self.checkfunctions]
        return "\n".join([chk for chk in checks if chk])

##########################
#     Syntax checks      #
##########################

    def checkfirstline(self):
        results = []
        firstlinedata = self.firstline.split()
        try:
            Nparticles = int(firstlinedata[0])
            if Nparticles != len(list(self.particlecounter.elements())):
                results.append(str(len(list(self.particlecounter.elements()))) + " particles in the event, but " + str(Nparticles) + " recorded in the first line! " + str(self.linenumber))
        except ValueError:
            results.append("Number of particles is " + firstlinedata[0] + ", not an integer! " + str(self.linenumber))
        try:
            processid = int(firstlinedata[1])
            if processid not in self.processidlist:
                results.append("Invalid process id %i!" % processid)
        except ValueError:
            results.append("Process id is " + firstlinedata[1] + ", not an integer! " + str(self.linenumber))
        try:
            eventweight = float(firstlinedata[2])
            if not usefulstuff.isfinite(eventweight): raise ValueError
        except ValueError:
            results.append("Event weight is " + firstlinedata[2] + ", not a number! " + str(self.linenumber))
        try:
            scale = float(firstlinedata[3])
            if not usefulstuff.isfinite(scale): raise ValueError
        except ValueError:
            results.append("Scale is " + firstlinedata[3] + ", not a number! " + str(self.linenumber))
        try:
            alphaQED = float(firstlinedata[4])
            if not usefulstuff.isfinite(alphaQED): raise ValueError
        except ValueError:
            results.append("alphaQED is " + firstlinedata[4] + ", not a number! " + str(self.linenumber))
        try:
            alphaQCD = float(firstlinedata[5])
            if not usefulstuff.isfinite(alphaQCD): raise ValueError
        except ValueError:
            results.append("alphaQCD is " + firstlinedata[5] + ", not a number! " + str(self.linenumber))
        return "\n".join(results)

    def checkstatus(self):
        results = []
        for p in self.particlelist:
            if p.status() == -999:
                continue                 #already gave an error
            if p.status() not in [-1, 1, -2, 2, 3, -9]:
                results.append("Particle " + str(p) + " has unknown status " + str(p.status()) +  "! " + str(self.linenumber))

            if (p.status() == -1 or p.status() == -9) and any(p.mothers()):
                results.append("Particle " + str(p) + " with status " + str(p.status()) + " has mothers " + str(p.mothers()) + "! " + str(self.linenumber))
            if p.status() != -1 and p.status() != -9 and not any(p.mothers()):
                results.append("Particle " + str(p) + " with status " + str(p.status()) + " has no mothers! " + str(self.linenumber))

            if p.status() == 1 and p.kids():
                results.append("Particle " + str(p) + " with status 1 has kids " + str(p.kids()) + "! " + str(self.linenumber))
            if p.status() != 1 and not p.kids():
                results.append("Particle " + str(p) + " with status " + str(p.status()) + " has no kids! " + str(self.linenumber))

            if p.status() == -2 and p.lhemass() > 0:
                results.append("Particle " + str(p) + " with status -2 has m = " + str(p.lhemass()) + " > 0! " + str(self.linenumber))
            if p.status() != -2 and p.status() != 3 and p.lhemass() < 0:
                results.append("Particle " + str(p) + " with status " + str(p.status()) + " has m = " + str(p.lhemass()) + " < 0! " + str(self.linenumber))
        return "\n".join(results)

    def checkinvmass(self):
        results = []
        for p in self.particlelist:
            if abs(p.invmass() - p.lhemass()) >= config.invmasstolerance:
                results.append("Mass is inconsistent! " + str(self.linenumber) + "(" + str(p) + ")\n" +
                               "invariant mass = " + str(p.invmass()) + "\n" +
                               "LHE mass       = " + str(p.lhemass()))
        return "\n".join(results)

    def checkPDGmass(self):
        results = []
        for p in self.particlelist:
            if p in config.checkPDGmasslist:
                if abs(p.lhemass() - p.PDGmass()) >= config.PDGmasstolerance * p.PDGmass():
                    results.append("Mass is wrong! " + str(self.linenumber) + "(" + str(p) + ")\n" +
                                   "PDG mass = " + str(p.PDGmass()) + "\n" +
                                   "LHE mass = " + str(p.lhemass()))
        return "\n".join(results)

    def checkmomentum(self):
        results = []
        for v in self.vertices.values():
            if v.momentumin() != v.momentumout():
                results.append("no momentum conservation! " + str(self.linenumber) + "\n" +
                               "mom momentum  = " + str(v.momentumin()) + str(v.particlesin()) + "\n" +
                               "kids momentum = " + str(v.momentumout()) + str(v.particlesout()))
        return "\n".join(results)

    def checkcharge(self):
        results = []
        for v in self.vertices.values():
            if v.chargein() != v.chargeout():
                results.append("no charge conservation! " + str(self.linenumber) + "\n" +
                               "mom charge  = " + str(v.chargein()) + str(v.particlesin()) + "\n" +
                               "kids charge = " + str(v.chargeout()) + str(v.particlesout()))
        return "\n".join(results)

    def checkcolor(self):
        results = []
        for p in self.particlelist:
            if not p.color() and ((p in globalvariables.globalvariables.quarks and p.id() > 0) or p in globalvariables.globalvariables.gluon):
                results.append(str(p) + " has no color! " + str(self.linenumber))
            if p.color() and not ((p in globalvariables.globalvariables.quarks and p.id() > 0) or p in globalvariables.globalvariables.gluon):
                results.append(str(p) + " has color " + str(p.color()) + "! " + str(self.linenumber))

            if not p.anticolor() and ((p in globalvariables.globalvariables.quarks and p.id() < 0) or p in globalvariables.globalvariables.gluon):
                results.append(str(p) + " has no anticolor! " + str(self.linenumber))
            if p.anticolor() and not ((p in globalvariables.globalvariables.quarks and p.id() < 0) or p in globalvariables.globalvariables.gluon):
                results.append(str(p) + " has anticolor " + str(p.anticolor()) + "! " + str(self.linenumber))

        for c in self.colors.values():
            if not c.check():
                results.append("color line " + str(c) + " doesn't make sense! " + str(self.linenumber) + "\n" +
                               "particles involved: " + str(c.particles.union(c.antiparticles)))
        return "\n".join(results)

    def is4e(self):
        return self.count(globalvariables.globalvariables.electrons) == 4
    def is4mu(self):
        return self.count(globalvariables.globalvariables.muons) == 4
    def is2e2mu(self):
        return self.count(globalvariables.globalvariables.electrons) == self.count(globalvariables.globalvariables.muons) == 2
