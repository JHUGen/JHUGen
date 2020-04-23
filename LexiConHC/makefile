CMSSWINCLUDEPATHS    = $(CMSSW_SEARCH_PATH)
CMSSWINCLUDES        = $(subst :, -I,$(CMSSWINCLUDEPATHS))
CMSSWFWLITEPATHS     = $(CMSSW_FWLITE_INCLUDE_PATH)
CMSSWFWLITEINCLUDES  = $(subst :, -I,$(CMSSWFWLITEPATHS))
BASEINCLUDE          = $(CMSSW_BASE)/src
MELADIR              = $(BASEINCLUDE)/ZZMatrixElement/MELA
CMSDATATOOLSDIR      = $(BASEINCLUDE)/CMSDataTools/AnalysisTree
CANDLORECASTERDIR    = $(BASEINCLUDE)/MelaAnalytics/CandidateLOCaster
EVENTCONTAINERDIR    = $(BASEINCLUDE)/MelaAnalytics/EventContainer
GENERICMECOMPUTERDIR = $(BASEINCLUDE)/MelaAnalytics/GenericMEComputer
COMPILEPATH          = $(PWD)/
INCLUDEDIR           = $(COMPILEPATH)interface/
SRCDIR               = $(COMPILEPATH)src/
BINDIR               = $(COMPILEPATH)bin/
OUTLIB               = $(COMPILEPATH)lib/
RUNDIR               = $(COMPILEPATH)

ROOTCFLAGS    = $(shell root-config --cflags) -Lrootlib -lEG
ROOTGLIBS     = $(shell root-config --glibs) -Lrootlib -lEG

CXX           = g++
CXXFLAGS      = -fPIC -g -O2

ARCH         := $(shell root-config --arch)

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
CXX           = c++ -lm
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O2 -j 4
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
UNDEFOPT      = suppress
endif

CXXFLAGS      += $(ROOTCFLAGS) -I$(ROOFITSYS)/include/ -I$(CMSSWINCLUDES) -I$(CMSSWFWLITEINCLUDES) -I$(INCLUDEDIR) -I$(MELADIR)/interface/ -I$(CMSDATATOOLSDIR)/interface/ -I$(CANDLORECASTERDIR)/interface/ -I$(EVENTCONTAINERDIR)/interface/ -I$(GENERICMECOMPUTERDIR)/interface/

NGLIBS         = $(ROOTGLIBS)
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L $(ROOFITSYS)/lib/ -lMinuit -lRooFitCore -lRooFit #-lboost_regex
# Libraries for common user packages
NGLIBS        += -L $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lZZMatrixElementMELA -lCMSDataToolsAnalysisTree -lMelaAnalyticsCandidateLOCaster -lMelaAnalyticsEventContainer -lMelaAnalyticsGenericMEComputer
# Libraries needed for CMSSW dependencies (mainly due to CMSDataTools)
NGLIBS        += -L $(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -lDataFormatsCommon -lDataFormatsMath
# Filter out libNew because it leads to floating-point exceptions in versions of ROOT prior to 6.08.02
# See the thread https://root-forum.cern.ch/t/linking-new-library-leads-to-a-floating-point-exception-at-startup/22404
GLIBS          = $(filter-out -lNew, $(NGLIBS))


SOURCES = $(wildcard $(SRCDIR)*.cc)
OBJECTS = $(SOURCES:.cc=.o)

.PHONY: all help compile clean

all: LHEAnalyzer trimPythia


# General rule for making object files
%.d:	$(SRCDIR)%.cc $(BINDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CXX) -MM -MT $@ -MT ${@:.d=.o} $(CXXFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(SRCDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cc $(BINDIR)%.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) $< -c -o $@ $(GLIBS)


LHEAnalyzer: $(OBJECTS) $(BINDIR)LHEAnalyzer.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)LHEAnalyzer $(OBJECTS) $(BINDIR)LHEAnalyzer.cc $(GLIBS)

trimPythia: $(OBJECTS) $(BINDIR)trimPythia.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)trimPythia $(OBJECTS) $(BINDIR)trimPythia.cc $(GLIBS) -lSimDataFormatsGeneratorProducts -lDataFormatsHepMCCandidate -lDataFormatsJetReco

clean:
	rm -rf $(OUTLIB)
	rm -rf $(RUNDIR)Pdfdata
	rm -f $(RUNDIR)input.DAT
	rm -f $(RUNDIR)process.DAT
	rm -f $(RUNDIR)ffwarn.dat
	rm -f $(RUNDIR)br.sm*
	rm -f $(SRCDIR)*.o
	rm -f $(SRCDIR)*.so
	rm -f $(SRCDIR)*.d
	rm -f $(BINDIR)*.o
	rm -f $(BINDIR)*.so
	rm -f $(BINDIR)*.d
	rm -f $(RUNDIR)*.cc
	rm -f $(RUNDIR)analyze*
	rm -f $(RUNDIR)LHEAnalyze*
	rm -f $(RUNDIR)trimPythia*
	rm -f $(RUNDIR)*.o
	rm -f $(RUNDIR)*.so
	rm -f $(RUNDIR)*.d
	rm -f $(RUNDIR)*.pcm
	rm -f $(RUNDIR)*.py
	rm -f $(RUNDIR)*.pyc
