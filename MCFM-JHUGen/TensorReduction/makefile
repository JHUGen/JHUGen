FTNCHEKPATH = /home/ellis/Fortran/Ftnchek/ftnchek-3.1.2
FORCHKPATH = /usr/local/bin/
FFLAGS 	= -g -Wall 
HERE = $(PWD)
LIBNAME = QCDLoop
LIBDIR = $(HERE)/Lib
INCDIR = $(HERE)/pv/Include
PVDIR  = $(HERE)/pv
OVDIR  = $(HERE)/ov
QLDIR  = $(HERE)/../QCDLoop/ql
FFDIR  = $(HERE)/../QCDLoop/ff
RRDIR  = $(HERE)/recur
ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif

# -----------------------------------------------------------------------------
# This makes the libraries that we need

libs:
	(cd $(RRDIR)/smallG; make) 
	(cd $(RRDIR)/smallY; make) 
	(cd $(RRDIR)/smallP; make) 
	(cd $(RRDIR)/smallF; make) 
	(cd $(PVDIR); make) 
	(cd $(OVDIR); make) 

clean:
	(cd $(RRDIR)/smallG; make clean) 
	(cd $(RRDIR)/smallY; make clean) 
	(cd $(RRDIR)/smallP; make clean) 
	(cd $(RRDIR)/smallF; make clean) 
	(cd $(PVDIR); make clean) 
	(cd $(OVDIR); make clean) 

# -----------------------------------------------------------------------------

# DO NOT DELETE


