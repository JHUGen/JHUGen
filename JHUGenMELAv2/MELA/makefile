# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else
ifndef ROOFITSYS
all:
	@echo "ROOTFITSYS is not set. Please set ROOT environment properly"; echo
else

RM  = /bin/rm

MELADIR = .
MELASRCDIR = $(MELADIR)/src
MELADATADIR = $(MELADIR)/data
# Modify MELALIBDIR for the gcc version as needed
MELALIBDIR = $(MELADATADIR)/slc6_amd64_gcc530

ROOTCFLAGS = $(shell root-config --cflags)
ROOTGLIBS = $(shell root-config --evelibs) -Lrootlib
NLIBS = $(ROOTGLIBS)
NLIBS += -L$(ROOFITSYS)/lib -lMinuit -lRooFitCore -lRooFit -lboost_regex -lboost_filesystem
LIBS = $(filter-out -lNew, $(NLIBS))

ROOFITINC = -I$(ROOFITSYS)/include
ROOTINC = $(ROOFITINC)

CC = g++
CPPINC = -I$(MELADIR)/interface $(ROOTINC)
CPPLOAD = -L$(MELALIBDIR) -lmcfm_703 -ljhugenmela $(LIBS)
CPPOPTS =  -fPIC -O2 -ftree-vectorize -fipa-pta -felide-constructors -fvisibility-inlines-hidden -fno-math-errno \
	--param vect-max-version-for-alias-checks=50 -fmessage-length=0 -fdiagnostics-show-option \
	 $(shell root-config --cflags)
CPPFLAGS = $(CPPOPTS) $(CPPINC) $(CPPLOAD)

LINKER = g++
LINKERFLAGS = -Wl,-rpath=$(MELALIBDIR),-soname,$(LIB) $(CPPLOAD)

SOURCES = $(wildcard $(MELASRCDIR)/*.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIB = libMELA.so

.PHONY: all help compile clean

all: SRCFILES LinkDefOut LinkDef $(LIB)

LinkDefOut: $(MELASRCDIR)/LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
        rootcint -f $(MELASRCDIR)/LinkDef_out.cc -c -p $(CPPINC) $(MELASRCDIR)/LinkDef.h; \
        cat $(MELASRCDIR)/LinkDef.h $(MELASRCDIR)/LinkDef_out.cc > $(MELASRCDIR)/LinkDef_out.cxx; \
	rm $(MELASRCDIR)/LinkDef_out.cc; mv $(MELASRCDIR)/LinkDef_out_rdict.pcm $(MELALIBDIR)/LinkDef_out_rdict.pcm

# General rule for making object files
%.d:	$(MELASRCDIR)/%.cc $(MELASRCDIR)/%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CPPFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(MELASRCDIR)/%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CPPFLAGS) $< -c -o $@

SRCFILES: $(OBJECTS)
	
LinkDef: $(MELASRCDIR)/LinkDef_out.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CPPOPTS) $(CPPINC) $< -c -o $(MELASRCDIR)/LinkDef_out.o

$(LIB):	$(OBJECTS)
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) $(MELASRCDIR)/LinkDef_out.o -o $@; \
	mv $(LIB) $(MELALIBDIR)/$(LIB)

clean:
	$(QUIET) rm -v -f \
	$(MELASRCDIR)/*.o $(MELASRCDIR)/*.d $(MELASRCDIR)/LinkDef_out* $(MELALIBDIR)/$(LIB) $(MELALIBDIR)/LinkDef_out_rdict.pcm; \
	echo "Done cleaning"


endif
endif

