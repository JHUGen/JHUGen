from fractions import Fraction
import particledata

class ParticleType(object):
    def __init__(self, particleorid):
        self.__id = int(particleorid)
        if particleorid < 0 and -particleorid in globalvariables.globalvariables.neutralbosons.ids():
            self.__id = -self.__id

    def id(self):
        return self.__id
    def __str__(self):
        if self.id() < 0:
            result = str(-self)
            if "-" in result:
                result = result.replace("-","+")
            elif "+" in result:
                result = result.replace("+","-")
            elif self in globalvariables.globalvariables.neutrinos:
                result = result.replace("nu","nubar")
            elif self in globalvariables.globalvariables.quarks:
                result = result + "bar"
            else:
                raise ValueError("Invalid particle id! " + str(self.id()))
            return result
        try:
            return particledata.particlename[self.id()]
        except KeyError:
            raise ValueError("Invalid particle id! " + str(self.id()))

    def __neg__(self):
        return ParticleType(-self.id())

    def charge(self):
        if self.id() < 0:
            return -(-self).charge()
        if self in globalvariables.globalvariables.downtypequarks:
            return Fraction(-1, 3)
        if self in globalvariables.globalvariables.uptypequarks:
            return Fraction(2, 3)
        if self in globalvariables.globalvariables.leptons:
            return -1
        if self in globalvariables.globalvariables.neutrinos:
            return 0
        if self in globalvariables.globalvariables.neutralbosons:
            return 0
        if self in globalvariables.globalvariables.W:
            return 1
        raise ValueError("Invalid particle id: " + str(self.id()))

    def PDGmass(self):
        if self.id() < 0:
            return (-self).PDGmass()
        return particledata.PDGmass[self.id()]
    def PDGmasserror(self):
        if self.id() < 0:
            return (-self).PDGmasserror()
        return particledata.PDGmasserror[self.id()]

    def baryonnumber(self):
        if self.id() < 0:
            return -(-self).baryonnumber()
        if self in globalvariables.globalvariables.quarks:
            return Fraction(1,3)
        return 0
    def leptonnumber(self, generation):
        if self.id() < 0:
            return -(-self).leptonnumber(generation)
        if (generation == 1 and self.id() in [11, 12]
         or generation == 2 and self.id() in [13, 14]
         or generation == 3 and self.id() in [15, 16]):
            return 1
        return 0

    def __hash__(self):
        return self.id()
    def __eq__(self, other):
        try:
            return self.id() == other.id()
        except AttributeError:
            return False
    def __ne__(self, other):
        return not self == other
    def __int__(self):
        return self.id()
    def __bool__(self):
        return True

import globalvariables
