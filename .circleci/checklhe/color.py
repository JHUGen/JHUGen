import usefulstuff

class Color(object):
    def __init__(self, id):
        self.id = id
        self.particles = usefulstuff.printableset()
        self.antiparticles = usefulstuff.printableset()

    def addparticle(self, p):
        if p.color() == self.id:
            self.particles.add(p)
        if p.anticolor() == self.id:
            self.antiparticles.add(p)

    def check(self):
        return linemakessense(self.particles, self.antiparticles)

    def __str__(self):
        return str(self.id)


def linemakessense(particles, antiparticles, start = None, end = None):
    particles = particles.copy()
    antiparticles = antiparticles.copy()
    if not particles.isdisjoint(antiparticles):
        return False

    if len(particles) == len(antiparticles) == 0:
        return True
    if start is not None and end is not None and start is end and len(particles) + len(antiparticles) == 1 and (start in particles or start in antiparticles):
        return True

    for p in particles:
        if p.startvertex is None:
            if start is not None and start is not p:
                return False
            start = p
        if p.endvertex is None:
            if end is not None and end is not p:
                return False
            end = p
    for p in antiparticles:
        if p.endvertex is None:
            if start is not None and start is not p:
                return False
            start = p
        if p.startvertex is None:
            if end is not None and end is not p:
                return False
            end = p
    if start is None and end is not None or start is not None and end is None:
        return False

    if start is None and end is None:
        try:
            start = list(particles)[0]
            end = list(particles)[0]
        except IndexError:     #no particles, only antiparticles
            start = list(antiparticles)[0]
            end = list(antiparticles)[0]

    if start in particles:
        nextvertex = start.endvertex
        if start is not end:
            particles.remove(start)
    elif start in antiparticles:
        nextvertex = start.startvertex
        if start is not end:
            antiparticles.remove(start)
    else:
        assert(0)

    if len(particles) == len(antiparticles) == 0:
        return True

    possiblenextparticles = usefulstuff.printableset()
    for p in particles:
        if p.startvertex is nextvertex:
            possiblenextparticles.add(p)
    for p in antiparticles:
        if p.endvertex is nextvertex:
            possiblenextparticles.add(p)

    return any(linemakessense(particles, antiparticles, nextparticle, end) for nextparticle in possiblenextparticles)
        

class Colors(dict):
    def __missing__(self, id):
        self[id] = Color(id)
        return self[id]
