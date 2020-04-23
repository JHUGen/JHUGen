COMPILEPATH          = $(PWD)/
INCLUDEDIR           = $(COMPILEPATH)interface/
SRCDIR               = $(COMPILEPATH)src/
BINDIR               = $(COMPILEPATH)bin/
OUTLIB               = $(COMPILEPATH)lib/
RUNDIR               = $(COMPILEPATH)

CXX           = g++
CXXFLAGS      = -fPIC -g -O2

CXXFLAGS      += -I$(INCLUDEDIR)

SOURCES = $(wildcard $(SRCDIR)*.cc)
OBJECTS = $(SOURCES:.cc=.o)

.PHONY: all help compile clean

all: LexiConHC


# General rule for making object files
%.d:	$(SRCDIR)%.cc $(BINDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CXX) -MM -MT $@ -MT ${@:.d=.o} $(CXXFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(SRCDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cc $(BINDIR)%.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) $< -c -o $@


LexiConHC: $(OBJECTS) $(BINDIR)LexiConHC.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)LexiConHC $(OBJECTS) $(BINDIR)LexiConHC.cc $(GLIBS)

clean:
	rm -rf $(OUTLIB)
	rm -f $(SRCDIR)*.o
	rm -f $(SRCDIR)*.so
	rm -f $(SRCDIR)*.d
	rm -f $(BINDIR)*.o
	rm -f $(BINDIR)*.so
	rm -f $(BINDIR)*.d
	rm -f $(RUNDIR)*.cc
	rm -f $(RUNDIR)LexiCon*
	rm -f $(RUNDIR)*.o
	rm -f $(RUNDIR)*.so
	rm -f $(RUNDIR)*.d
	rm -f $(RUNDIR)*.pcm
	rm -f $(RUNDIR)*.py
	rm -f $(RUNDIR)*.pyc
