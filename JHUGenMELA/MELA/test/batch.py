#!/usr/bin/env python

import os
import sys
import commands

from ROOT import TFile

def processDirectory ( args, dirname, filenames ):
    print "processing " + dirname
    for filename in filenames:
        fullname = dirname + "/" + filename
        
        flavor = 10
        if ("4e" in dirname):
            flavor = 1
        if ("4mu" in dirname):
            flavor = 2
        if ("2mu2e" in dirname) or ("2e2mu" in dirname):
            flavor = 3

        if ("data" in  dirname) and ("DoubleEle" in filename):
            flavor = 1
        if ("data" in  dirname) and ("DoubleMu" in filename):
            flavor = 2
        if ("data" in  dirname) and ("DoubleOr" in filename):
            flavor = 3

        if ("CR" in  dirname) and ("DoubleEle" in filename):
            flavor = 1
        if ("CR" in  dirname) and ("DoubleMu" in filename):
            flavor = 2
        if ("CR" in  dirname) and ("DoubleOr" in filename):
            flavor = 3

        if ("CR" in dirname) and ("EEEE" in filename):
            flavor = 1
        if ("CR" in dirname) and ("MMMM" in filename):
            flavor = 2
        if ("CR" in dirname) and ("EEMM" in filename):
            flavor = 3
        if ("CR" in dirname) and ("MMEE" in filename):
            flavor = 3
                            

        if ("withProbabilities" in filename):
            flavor = 10  # already processed
                        
        if not(".root" in filename):
            flavor = 10  # only process root files
            
        sqrts=7
        if ("8TeV" in dirname):
            sqrts = 8
            


        print " " * 4 + filename + " with flavor " + str(flavor) + " and sqrts = " + str(sqrts) 

        

        if flavor!=10: # looks like a valid file, prepare string
            command = "root -q -b  addProbtoTree.C\\(\\\"" + fullname[:-5] + "\\\","+str(flavor)+",-1,"+str(sqrts)+"\\)\n"
            #create batch script
            commands.getstatusoutput("cp batchscript.csh batchscript_tmp.csh")
            file = open('batchscript_tmp.csh', 'a')
            file.write(command)
            file.close()
            commands.getstatusoutput("bsub -q 8nh < batchscript_tmp.csh" )
            #exit(0)
            

def main():
    base_dir = "."
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    file = open('./batchscript.csh', 'w')
    file.write("cd " + sys.argv[2] + "\n")
    file.write("pwd\n")
    file.write("eval `scramv1 runtime -sh`\n")
    file.write("cd - \n")
    file.write("cp " + sys.argv[2] + "/addProbtoTree.C .  \n")
    file.close()
    
    os.path.walk( base_dir, processDirectory, None )

if __name__ == "__main__":
    main()
