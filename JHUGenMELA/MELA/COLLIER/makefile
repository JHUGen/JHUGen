MAINDIR = .
AUXDIR = $(MAINDIR)/Aux
COLIDIR = $(MAINDIR)/COLI
TENSORSDIR = $(MAINDIR)/tensors
DDLIBDIR = $(MAINDIR)/DDlib
#MODULES=$(MAINDIR)/modules
#OBJECTS=$(MAINDIR)/objects
MODULES=$(MAINDIR)
OBJECTS=$(MAINDIR)
ROOFITINCLUDE = 
RM  = /bin/rm
INCLUDE = -I$(MAINDIR) $(ROOFITINCLUDE)
Comp = gfort

ifeq ($(Comp),ifort)
	fcomp = ifort -132 -O2-lifcore -fPIC
endif
ifeq ($(Comp),gfort)
	fcomp = gfortran  -ffixed-line-length-132  -fno-default-integer-8 -fPIC -O2 -funroll-loops -Wtabs
endif

ifeq ($(Comp),gfort)
	ccomp = gcc
endif
ifeq ($(Comp),ifort)
	ccomp = icc
endif

ifeq ($(Comp),gfort)
	clib =  -lm -lgfortran
endif
ifeq ($(Comp),ifort)
	clib =  -lm -lirc
endif


LIB = libcollier.so

Deps = $(MAINDIR)/collier_global.F90 \
	$(AUXDIR)/Combinatorics.F90 \
	$(AUXDIR)/master.F90 \
	$(AUXDIR)/cache.F90 \
	$(COLIDIR)/coli_aux2.F90 \
	$(COLIDIR)/coli_aux.F \
	$(COLIDIR)/coli_stat.F90 \
	$(MAINDIR)/collier_aux.F90 \
	$(TENSORSDIR)/InitTensors.F90 \
	$(MAINDIR)/collier_init.F90 \
	$(COLIDIR)/reductionAB.F90 \
	$(COLIDIR)/reductionC.F90 \
	$(COLIDIR)/reductionD.F90 \
	$(COLIDIR)/reductionEFG.F90 \
	$(COLIDIR)/reductionTN.F90 \
	$(MAINDIR)/collier_coefs.F90 \
	$(TENSORSDIR)/BuildTensors.F90 \
	$(TENSORSDIR)/TensorReduction.F90 \
	$(MAINDIR)/collier_tensors.F90 \
	$(DDLIBDIR)/DD_global.F90 \
	$(MAINDIR)/COLLIER.F90 \
	$(COLIDIR)/coli_b0.F \
	$(COLIDIR)/coli_d0.F \
	$(COLIDIR)/coli_d0reg.F \
	$(COLIDIR)/coli_c0.F \
	$(DDLIBDIR)/DD_aux.F \
	$(DDLIBDIR)/DD_2pt.F \
	$(DDLIBDIR)/DD_3pt.F \
	$(DDLIBDIR)/DD_4pt.F \
	$(DDLIBDIR)/DD_5pt.F \
	$(DDLIBDIR)/DD_6pt.F \
	$(DDLIBDIR)/DD_to_COLLIER.F

Objs = collier_global.o \
	Aux_Combinatorics.o \
	Aux_master.o \
	Aux_cache.o \
	COLI_coli_aux2.o \
	COLI_coli_aux.o \
	COLI_coli_stat.o \
	collier_aux.o \
	tensors_InitTensors.o \
	collier_init.o \
	COLI_reductionAB.o \
	COLI_reductionC.o \
	COLI_reductionD.o \
	COLI_reductionEFG.o \
	COLI_reductionTN.o \
	collier_coefs.o \
	tensors_BuildTensors.o \
	tensors_TensorReduction.o \
	collier_tensors.o \
	DDlib_DD_global.o \
	COLLIER.o \
	COLI_coli_b0.o \
	COLI_coli_d0.o \
	COLI_coli_d0reg.o \
	COLI_coli_c0.o \
	DDlib_DD_aux.o \
	DDlib_DD_2pt.o \
	DDlib_DD_3pt.o \
	DDlib_DD_4pt.o \
	DDlib_DD_5pt.o \
	DDlib_DD_6pt.o \
	DDlib_DD_to_COLLIER.o


${Objs}: $(Deps)
	@echo " "
	@echo " Compiling COLLIER dependencies"
	$(fcomp) -c $(Deps)
	@echo " "
	@echo " Compiling COLLIER library"
	g++ -Wl,-soname,$(LIB) -shared -o $(LIB) *.o

clean:
	@echo " deleting object files"
	rm -f *.so $(OBJECTS)/*.o $(MODULES)/*.mod


# supresses command calls
.SILENT:
