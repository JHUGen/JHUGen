ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
HERE   = $(PWD)
INCDIR  = $(HERE)/Include/
FFLAGS =  -fno-f2c -fopenmp -O2 -g -I $(INCDIR)
QLDIR  = $(HERE)/ql
FFDIR  = $(HERE)/ff
QPDIR  = $(HERE)/ql/qp/
QFDIR  = $(HERE)/ff/qf/

# -----------------------------------------------------------------------------
# Specify the object files. 

TEST= \
test.o 

QPTEST= \
qptest.o 



.f.o:	
		$(FC) $(FFLAGS) -c $*.f


# -----------------------------------------------------------------------------

all: test

qptest: qpdir qfdir qldir ffdir		\
	$(QPTEST)
	$(FC) $(FFLAGS) -L$(QPDIR) -L$(QFDIR) -L$(QLDIR) -L$(FFDIR) -o $@ \
            $(QPTEST) -lqpqcdloop -lqf -lqcdloop -lff 

test: qldir ffdir 		\
	$(TEST)
	$(FC) $(FFLAGS) -L$(QLDIR) -L$(FFDIR) -o $@ \
            $(TEST) -lqcdloop -lff 

ffdir:
	(cd $(FFDIR); make; ranlib libff.a) 

qfdir:
	(cd $(QFDIR); make; ranlib libqf.a) 

qpdir:
	(cd $(QPDIR); make; ranlib libqpqcdloop.a) 

qldir:
	(cd $(QLDIR); make; ranlib libqcdloop.a) 


# -----------------------------------------------------------------------------

clean:
	- rm -f *.o *.s *.prj *~ core 

# -----------------------------------------------------------------------------

# DO NOT DELETE


