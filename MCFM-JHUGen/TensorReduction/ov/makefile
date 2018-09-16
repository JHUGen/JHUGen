# Makefile routine.
VPATH		= $(PWD)/obj
INCPATH  	= $(PWD)/../Include
OUTPUT_OPTION	= -o $(PWD)/obj/$@

OLODIR          = /Users/johnmc/Research/OneLOop-3.3.1

ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
# For linking against OneLOop
FFLAGS 	= -O2 -fno-f2c -g -I$(INCPATH) -Iobj -Jobj -I$(OLODIR)
# no linking
FFLAGS 	= -O2 -fno-f2c -g -I$(INCPATH) -Iobj -Jobj


LIBRARY = libov.a

# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
trI1.o \
trI2.o \
trI3.o \
trI4.o \
TRsettensorcontrol.o \
del2.o \
del3.o \
del4.o \
doBtensor.o \
doCtensor.o \
doDtensor.o \
doEtensor.o \
ludcmp.o \
mprove.o \
ovArraysetup.o \
ovBtensor.o \
ovCcheck.o \
ovCtensor.o \
ovDcheck.o \
ovDtensor.o \
ovE0scalar.o \
ovEtensor.o \
ovEcheck.o \
ovIpiksrt.o \
ovdot.o \
ovdsvdcmp.o \
ovdpythag.o \
ovswitch.o \
ovw2.o \
ovw3.o \
pvfndd.o \
zlubksb.o \
zsvbksb.o
      

      
$(LIBRARY):	$(OBJS)
		ar cru $(LIBRARY) $(patsubst %,obj/%,$(OBJS))
		ranlib $(LIBRARY)
 
# -----------------------------------------------------------------------------
# Specify other options.


clean:
	- rm -f *.o obj/*.o obj_omp/*.o *.s *.prj *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE


# LocalWords:  Makefile
