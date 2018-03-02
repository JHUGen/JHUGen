#! /usr/bin/env python
import sys

#https://root.cern.ch/phpBB3/viewtopic.php?t=3198
tempargv = sys.argv
sys.argv.insert(0, '-b')
import ROOT
del sys.argv[0]
#Try both ways
ROOT.gROOT.SetBatch(True)

import lhefile

totalevents = 0
total4e = 0
total4mu = 0
total2e2mu = 0

for file in sys.argv[1:]:
    print file
    with lhefile.LHEFile(file) as f:
        for ev in f:
            check = ev.check()
            if check:
                f.raiseerror(check)
        totalevents += f.nevents
        total4e += f.n4e
        total4mu += f.n4mu
        total2e2mu += f.n2e2mu

if len(sys.argv) >= 3: #at least 2 files
    print totalevents, "total events"
    print total4e, "total 4e events"
    print total4mu, "total 4mu events"
    print total2e2mu, "total 2e2mu events"
