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


MELADIR = .
MELASRCDIR = $(MELADIR)/src
MELADATADIR = $(MELADIR)/data
MELASLCDIR = $(MELADATADIR)/slc6_amd64_gcc530
ROOFITINCLUDE = 
RM  = /bin/rm
CPPINC = -I$(MELADIR)/interface $(ROOFITINCLUDE)

CC = g++
CPPFLAGS = -fPIC -ansi -D_GNU_SOURCE -g -O2  $(shell root-config --cflags) $(CPPINC)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
ifdef CMSSW_RELEASE_BASE
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -L$(CMSSW_RELEASE_BASE)/external/$(SCRAM_ARCH)/lib $(shell root-config --libs) -lEG -lGenVector
else
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif
endif


SOURCES = $(wildcard $(MELASRCDIR)/*.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIB = libMELA.so

.PHONY: all help compile clean

all: SRCFILES LinkDefOut LinkDef $(LIB)

LinkDefOut: $(MELASRCDIR)/LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
        rootcint -f $(MELASRCDIR)/LinkDef_out.cc -c -p $(CPPINC) $(MELASRCDIR)/LinkDef.h; \
        cat $(MELASRCDIR)/LinkDef.h $(MELASRCDIR)/LinkDef_out.cc > $(MELASRCDIR)/LinkDef_out.cxx; rm $(MELASRCDIR)/LinkDef_out.cc

# General rule for making object files
%.d:	$(MELASRCDIR)/%.cc $(MELASRCDIR)/%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CPPFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(MELASRCDIR)/%.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CPPFLAGS) $< -c -o $@

%.o: 	$(MELASRCDIR)/%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CPPFLAGS) $< -c -o $@

SRCFILES: $(OBJECTS)
	
LinkDef: $(MELASRCDIR)/LinkDef_out.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CPPFLAGS) $< -c -o $(MELASRCDIR)/LinkDef_out.o

$(LIB):	$(OBJECTS)
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) $(MELASRCDIR)/LinkDef_out.o -o $@; \
	mv $(LIB) $(MELASLCDIR)/$(LIB)

clean:   
	$(QUIET) rm -v -f \
	$(MELASRCDIR)/*.o $(MELASRCDIR)/*.d $(MELASRCDIR)/LinkDef_out* $(MELASLCDIR)/$(LIB); echo "Done"


endif
endif

