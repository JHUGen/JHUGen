import globalvariables
import ROOT

#syntax checks
checkfirstline = True         #check that the first line in each event has the right syntax, and that the number of particles is correct
checkstatus = True            #check that the particle status is correct

#kinematic checks
checkinvmass = True           #check that E^2-\vec{p}^2=m^2, using all the numbers from the LHE
checkPDGmass = False          #check that the LHE mass = the PDG mass for the particle types in checkPDGmasslist
                              #   (configured below)
checkmomentum = True          #check momentum conservation in particle decays
checkcharge = True            #check charge conservation in particle decays
checkcolor = True             #check that the color lines make sense

#decay counts and checks
checkZZorWWassignment = True  #check that, in symmetric decays, the Zs/Ws are assigned so that one has mass as close as possible to m_Z/m_W

#Other options
raiseerror = False                #if this is true, raise an IOError when something is wrong
                                  # otherwise just print it
allowimplicithiggs = True         #Count ZZ or WW produced by the initial particles as a Higgs
                                  # Useful for MCFM background
ROOT.gErrorIgnoreLevel=ROOT.kError#Supress ROOT warnings, the main bad one is when eta = infinity


#tolerances
momentumtolerance = 1e-4      #GeV
invmasstolerance = 7e-2       #GeV
PDGmasstolerance = 0.05       #relative


#mostly internal, until the row of #####
startedinit = False
finishedinit = False

def init():
    global startedinit, finishedinit
    global checkPDGmasslist
    if not startedinit:
        startedinit = True
        globalvariables.init()
####################################################
#the stuff in here is configurable
        checkPDGmasslist = globalvariables.globalvariables.leptons
####################################################
        finishedinit = True

