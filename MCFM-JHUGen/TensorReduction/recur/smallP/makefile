# Makefile routine.
VPATH		= $(PWD)/obj
INCPATH1 = $(PWD)/../Include/
INCPATH2 = $(PWD)/../../Include/
OUTPUT_OPTION	= -o $(PWD)/obj/$@ 

ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
FFLAGS 	= -O2 -fno-f2c -g -I$(INCPATH1) -I$(INCPATH2) -Iobj -Jobj

LIBRARY	      = libsmallP.a


# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
Cfill_recur3.o \
Dfill_recur3.o \
runCP_0.o \
runCP_00.o \
runCP_0000.o \
runCP_0000i.o \
runCP_00i.o \
runCP_00ii.o \
runCP_00iii.o \
runCP_00iiii.o \
runCP_i.o \
runCP_ii.o \
runCP_iii.o \
runCP_iiii.o \
runCP_iiiii.o \
runP_0.o \
runP_00.o \
runP_0000.o \
runP_0000i.o \
runP_00i.o \
runP_00ii.o \
runP_00iii.o \
runP_00iiii.o \
runP_i.o \
runP_ii.o \
runP_iii.o \
runP_iiii.o \
runP_iiiii.o


$(LIBRARY):	$(OBJS) 
		ar cru $(LIBRARY) $(patsubst %,obj/%,$(OBJS))
		ranlib $(LIBRARY)
		mv $(LIBRARY) ..
 

# -----------------------------------------------------------------------------



all:		$(LIBRARY)


# -----------------------------------------------------------------------------

clean:
	- rm -f obj/*.o *.s *.prj *.a *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE


