import collections
import ROOT
import config
import globalvariables
import event
import particle

class LHEFile:
    def __init__(self, filename):
        globalvariables.init()
        if filename[-4:] != ".lhe":
            raise ValueError(filename + " does not end in .lhe")
        self.filename = filename
        self.f = open(filename)
        self.nevents = 0
        self.n4e = 0
        self.n4mu = 0
        self.n2e2mu = 0
        self.linenumber = 0
        self.incomment = False
        self.sawinitblock = False
        self.processidlist = []
        self.VegasNc2 = None

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        print "   ", self.nevents, "events"
        print "   ", self.n4e, "4e events"
        print "   ", self.n4mu, "4mu events"
        print "   ", self.n2e2mu, "2e2mu events"
        self.f.close()
        if self.VegasNc2 is not None and self.nevents != self.VegasNc2:
            self.raiseerror("VegasNc2={}, but {} events!".format(self.VegasNc2, self.nevents))

    def raiseerror(self, msg):
        if config.raiseerror:
            raise IOError(msg)
        else:
            print msg

    def readevent(self):
        while "<event" not in self.nextline():
            if not self.line:     #eof
                return None

            if "</event>" in self.line:
                self.raiseerror("Extra </event>! " + str(self.linenumber))

            if "<init>" in self.line:
                if self.sawinitblock:
                    self.raiseerror("Extra init block! " + str(self.linenumber))
                self.sawinitblock = True
                data = self.nextline().split()
                try:
                    [int(data[i]) for i in (0, 1, 4, 5, 6, 7, 8, 9)]
                    [float(data[i]) for i in (2, 3)]
                except (ValueError, IndexError):
                    self.raiseerror("Bad init line 1!")
                nprocesses = int(data[9])
                for p in range(nprocesses):
                    data = self.nextline().split()
                    if "</init>" in self.line:
                        self.raiseerror("Not enough lines in init block!")
                        break
                    try:
                        [float(data[i]) for i in (0, 1, 2)]
                        if float(data[0]) < 0:
                            self.raiseerror("Pythia doesn't like negative cross sections!")
                        int(data[3])
                        for i in range(3, len(data)):
                            self.processidlist.append(int(data[i]))
                    except (ValueError, IndexError):
                        self.raiseerror("Bad init line %i!" % (2+p))
                while "</init>" not in self.nextline() and "<event>" not in self.line:
                    if self.line.split():
                        self.raiseerror("Extra line in init block!")
                if "<event>" in self.line:
                    self.raiseerror("No </init>!")
                    break

            if "<!--" in self.line:
                if not self.line.strip().startswith("<!--"):
                    self.raiseerror("Warning: comment begins in the middle of a line\n"
                                    "(ok in itself, but other problems may not be detected)! " + str(self.linenumber))
                if self.incomment:
                    self.raiseerror("<!-- inside a comment! " + str(self.linenumber))
                self.line = self.line.replace("<!--", "", 1)
                if "<!--" in self.line:
                    self.raiseerror("Warning: multiple <!-- in one line\n"
                                    "(ok in itself, but other problems may not be detected!" + str(self.linenumber))
                self.incomment = True
            if "-->" in self.line:
                if not self.line.strip().endswith("-->"):
                    self.raiseerror("Warning: comment ends in the middle of a line\n"
                                    "(ok in itself, but problems may not be detected)! " + str(self.linenumber))
                if not self.incomment:
                    self.raiseerror("--> not preceded by <!--! " + str(self.linenumber))
                self.line = self.line.replace("-->", "", 1)
                if "-->" in self.line:
                    self.raiseerror("Warning: multiple --> in one line\n"
                                    "(ok in itself, but other problems may not be detected!" + str(self.linenumber))
                self.incomment = False
            if "--" in self.line and self.incomment:
                self.raiseerror("-- in a comment! " + str(self.linenumber))

            if self.incomment and "VegasNc2=" in self.line and ("VBFoffsh_run=*" in self.line or not any("Process={}".format(_) in self.line for _ in (66,67,68,69))):
                for argument in self.line.split():
                    if argument.startswith("VegasNc2="):
                        self.VegasNc2 = int(argument.split("=")[-1])


        if not self.sawinitblock:
            self.raiseerror("No <init>!")
        ev = event.Event(self.linenumber, self.processidlist)
        ev.setfirstline(self.nextline())
        while "</event>" not in self.nextline():
            if not self.line:
                self.raiseerror("File ends in the middle of an event!")
                return None
            if "<event" in self.line:
                self.raiseerror("Extra <event>! " + str(self.linenumber))
            try:
                ev.addparticle(self.line)
            except particle.BadParticleLineError:
                continue
        ev.finished()
        self.nevents += 1
        if ev.is4e(): self.n4e += 1
        if ev.is4mu(): self.n4mu += 1
        if ev.is2e2mu(): self.n2e2mu += 1
        return ev

    def nextline(self):
        self.linenumber += 1
        self.line = self.f.readline()
        return self.line

    def __iter__(self):
        return self
    def next(self):
        ev = self.readevent()
        if ev is not None:
            return ev
        raise StopIteration
