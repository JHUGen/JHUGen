ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
	false
else

MAINDIR = ./
#MODULES=$(MAINDIR)modules/
#OBJECTS=$(MAINDIR)objects/
MODULES=$(MAINDIR)
OBJECTS=$(MAINDIR)
ROOFITINCLUDE = 
RM  = /bin/rm
INCLUDE = -I$(MAINDIR) $(ROOFITINCLUDE)
Comp = gfort

# Linking LHAPDF
# link pdfs via LHA library ('Yes' or 'No')
UseLHAPDF=No
# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
MyLHADir=/./
# remember to export
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}
ifeq ($(UseLHAPDF),Yes)
   LHAPDFflags = -L$(MyLHADir) -lLHAPDF -DuseLHAPDF=1
else
   LHAPDFflags = -DuseLHAPDF=0
endif

# Linking the Collier library
UseCOLLIER=Yes
# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
MyCOLLIERDir=../data/${SCRAM_ARCH}/
MyCOLLIERInc=../COLLIER/
# remember to export
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}
ifeq ($(UseCOLLIER),Yes)
   COLLIERflags = -L$(MyCOLLIERDir) -lcollier -DuseCollier=1 -I$(MyCOLLIERInc)
else
   COLLIERflags = -DuseCollier=0
endif


ifeq ($(Comp),ifort)
	fcomp = ifort -fpp -O2 -vec-report0 -Dcompiler=1 -lifcore $(LHAPDFflags) $(COLLIERflags) -fPIC
endif
ifeq ($(Comp),gfort)
	fcomp = gfortran -O0 -ffree-line-length-none -Dcompiler=2 $(LHAPDFflags) $(COLLIERflags) -fno-automatic -fno-f2c  -fPIC -g 
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


LIB = libjhugenmela.so

all:    libs
	g++ -Wl,-soname,$(LIB) -shared -o $(LIB) *.o
	@echo " "
	@echo " WARNING: if you updated the fortran code,"
	@echo "          you need to move libjhugenmela.so"
	@echo '          to the appropriate folder in data/'
	@echo "          and then scram b.  Otherwise your changes"
	@echo "          will not apply to any code you run."
	@echo "          If you are running setup.sh this is done automatically."
#	mv *.o $(OBJECTS)
#	mv *.mod $(MODULES)

libs: CALLING_cpHTO.o \
	Cteq61Pdf.o mstwpdf.o NNPDFDriver.o \
	mod_Parameters.o mod_Misc.o mod_HashCollection.o mod_Kinematics.o \
	mod_JHUGenMELA.o mod_JHUGen.o \
	mod_Higgs.o mod_Graviton.o mod_Zprime.o \
	mod_HiggsJJ.o mod_HiggsJ.o mod_VHiggs.o \
	mod_TopDecay.o mod_TTBHiggs.o mod_THiggs.o

CALLING_cpHTO.o: ./CPS/CALLING_cpHTO.f
	@echo " "
	@echo " compiling CALLING_cpHTO.f with "$(Comp)
	$(fcomp) -c ./CPS/CALLING_cpHTO.f


Cteq61Pdf.o: ./pdfs/Cteq61Pdf.f
	@echo " "
	@echo " compiling Cteq61Pdf.f with "$(Comp)
	$(fcomp) -c ./pdfs/Cteq61Pdf.f


mstwpdf.o: ./pdfs/mstwpdf.f
	@echo " "
	@echo " compiling mstwpdf.f with "$(Comp)
	$(fcomp) -c ./pdfs/mstwpdf.f


NNPDFDriver.o: ./pdfs/NNPDFDriver.f
	@echo " "
	@echo " compiling NNPDFDriver.f with "$(Comp)
	$(fcomp) -c ./pdfs/NNPDFDriver.f


vegas.o: ./vegas.f
	@echo " "
	@echo " compiling vegas.f with "$(Comp)
	$(fcomp) -c ./vegas.f


mod_Parameters.o: mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Parameters.F90 with "$(Comp)
	$(fcomp) -c mod_Parameters.F90

mod_HashCollection.o: mod_HashCollection.F90
	@echo " "
	@echo " compiling mod_HashCollection.F90 with "$(Comp)
	$(fcomp) -c mod_HashCollection.F90

mod_JHUGenMELA.o: mod_JHUGenMELA.F90
	@echo " "
	@echo " compiling mod_JHUGenMELA.F90 with "$(Comp)
	$(fcomp) -c mod_JHUGenMELA.F90

mod_Misc.o: mod_Misc.F90
	@echo " "
	@echo " compiling mod_Misc.F90 with "$(Comp)
	$(fcomp) -c mod_Misc.F90

mod_Kinematics.o: mod_Kinematics.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Kinematics.F90 with "$(Comp)
	$(fcomp) -c mod_Kinematics.F90

mod_JHUGen.o: mod_JHUGen.F90 mod_Kinematics.F90 mod_Parameters.F90 mod_Misc.F90
	@echo " "
	@echo " compiling mod_JHUGen.F90 with "$(Comp)
	$(fcomp) -c mod_JHUGen.F90


mod_Higgs.o: mod_Higgs.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Higgs.F90 with "$(Comp)
	$(fcomp) -c mod_Higgs.F90


mod_Zprime.o: mod_Zprime.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Zprime.F90 with "$(Comp)
	$(fcomp) -c mod_Zprime.F90

mod_Graviton.o: mod_Graviton.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_Graviton.F90 with "$(Comp)
	$(fcomp) -c mod_Graviton.F90

mod_HiggsJ.o: mod_HiggsJ.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_HiggsJ.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJ.F90

mod_HiggsJJ.o: mod_HiggsJJ.F90 mod_Misc.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_HiggsJJ.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJJ.F90

mod_VHiggs.o: mod_VHiggs.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_VHiggs.F90 with "$(Comp)
	$(fcomp) -c mod_VHiggs.F90 

mod_TopDecay.o: mod_TopDecay.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_TopDecay.F90 with "$(Comp)
	$(fcomp) -c mod_TopDecay.F90  

mod_TTBHiggs.o: mod_TTBHiggs.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_TTBHiggs.F90 with "$(Comp)
	$(fcomp) -c mod_TTBHiggs.F90  

mod_THiggs.o: mod_THiggs.F90 mod_Parameters.F90
	@echo " "
	@echo " compiling mod_THiggs.F90 with "$(Comp)
	$(fcomp) -c mod_THiggs.F90  

clean:
	@echo " deleting object files"
	rm -f *.so $(OBJECTS)*.o $(MODULES)*.mod

endif

# supresses command calls
.SILENT:
