import particletype

#this should be a frozenset but it doesn't work for some reason
class ParticleCategory(set):
    def __init__(self, lst, Csymmetric = True):
        newlist = []
        for p in lst:
            type = particletype.ParticleType(p)
            newlist.append(type)
            if Csymmetric:
                newlist.append(-type)
        super(ParticleCategory, self).__init__(newlist)
    def __contains__(self, particle):
        try:
            if super(ParticleCategory, self).__contains__(particle):
                return True
        except TypeError:
            pass
        try:
            return super(ParticleCategory, self).__contains__(particletype.ParticleType(particle))
        except ValueError:
            return False
    def __str__(self):
        return " ".join(str(p) for p in self)
    def ids(self):
        return [p.id() for p in self]
    def __hash__(self):
        return hash(frozenset(self))
