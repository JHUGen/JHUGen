from math import isinf, isnan

class printablelist(list):
    def __str__(self, joiner = ", "):
        return "[" + joiner.join(str(a) for a in self) + "]"
    def __getitem__(self, key):
        if isinstance(key, slice):
            return printablelist(super(printablelist, self).__getitem__(key))
        return super(printablelist,self).__getitem__(key)
    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))
    def __getattr__(self, name):
        return printablelist([(a if a is None else a.__getattribute__(name)) for a in self])
    def __call__(self, *args, **kwargs):
        return printablelist([(a if a is None else a(*args, **kwargs)) for a in self])

class printableset(set):
    def __str__(self, joiner = ", "):
        return "(" + joiner.join(str(a) for a in self) + ")"
    def __repr__(self, joiner = ", "):
        return "set([" + joiner.join(str(a) for a in self) + "])"
    def __getattr__(self, name):
        return printablelist([(a if a is None else a.__getattribute__(name)) for a in self])   #note that this is a list!
    def __call__(self, *args, **kwargs):
        return printablelist([(a if a is None else a(*args, **kwargs)) for a in self])         #note that this is a list!

class printablefrozenset(frozenset):
    def __str__(self, joiner = ", "):
        return "(" + joiner.join(str(a) for a in self) + ")"
    def __repr__(self, joiner = ", "):
        return "frozenset([" + joiner.join(str(a) for a in self) + "])"
    def __getattr__(self, name):
        return printablelist([(a if a is None else a.__getattribute__(name)) for a in self])   #note that this is a list!
    def __call__(self, *args, **kwargs):
        return printablelist([(a if a is None else a(*args, **kwargs)) for a in self])         #note that this is a list!

class printabledict(dict):
    def __str__(self, joiner = ", "):
        return "{" + joiner.join(str(a) + ": " + str(self[a]) for a in self) + "}"
    def __getattr__(self, name):
        return printabledict({a: (self[a] if self[a] is None else self[a].__getattribute__(name)) for a in self})
    def __call__(self, *args, **kwargs):
        return printabledict({a: (self[a] if self[a] is None else self[a](*args, **kwargs)) for a in self})

def isfinite(number):
    return not isinf(number) and not isnan(number)
