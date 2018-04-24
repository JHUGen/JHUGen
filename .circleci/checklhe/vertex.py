import usefulstuff
import momentum

class Vertex(object):
    def __init__(self, in1, in2):
        self.__in = usefulstuff.printablefrozenset([in1, in2])   #If it's a decay (not a collision), in1 is the same as in2 so it's only entered in the set once
        self.__out = usefulstuff.printableset()

    def __hash__(self):
        return hash(tuple(a for a in self.__in))

    def __str__(self):
        return str(self.__in) + " --> " + str(self.__out)

    def addkid(self, kid):
        if not set(self.__in for mom in kid.mothers()):
            raise ValueError("Trying to add particle '%s' to a vertex '%s', but its mothers '%s' are different from the vertex!" % kid, self.__in, kid.mothers())
        self.__out.add(kid)

    def particlesin(self):
        return self.__in
    def particlesout(self):
        return self.__out

    def momentumin(self):
        if all(a is None for a in self.__in):
            return self.momentumout()               #momentum is automatically "conserved"
        return sum(self.__in.momentum(), momentum.Momentum(None, 0, 0, 0, 0))
    def momentumout(self):
        return sum(self.__out.momentum(), momentum.Momentum(None, 0, 0, 0, 0))
    def chargein(self):
        if all(a is None for a in self.__in):
            return self.chargeout()                 #charge is automatically "conserved"
        return sum(self.__in.charge())
    def chargeout(self):
        return sum(self.__out.charge())

class Vertices(dict):
    def __missing__(self, item):
        if all(a is None for a in item):
            return None
        itemlist = list(item)
        if len(itemlist) == 1:
            itemlist.append(itemlist[0])
        if not len(itemlist) == 2:
            raise ValueError("The only items that should be added to Vertices are collections of 1 or 2 particles")
        self[item] = Vertex(itemlist[0], itemlist[1])
        return self[item]
