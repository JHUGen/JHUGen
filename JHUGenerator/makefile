Here = $(PWD)
Makefile=makefile
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas

# compiler options, Comp='ifort' or 'gfort'
Opt  = No
Comp = ifort

# link pdfs via LHA library ('Yes' or 'No')
UseLHAPDF=No

# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
MyLHADir=/home/schulze/lib/LHAPDF-6.1.5/lib/
# remember to export 
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}



ifeq ($(Comp),gfort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O3
 else
    OptOpts  = -O0
 endif
 ifeq ($(UseLHAPDF),Yes)
    LHAPDFflags = -L$(MyLHADir) -lLHAPDF -DuseLHAPDF=1
 else
    LHAPDFflags = -DuseLHAPDF=0
 endif 
 fcomp = f95 $(OptOpts) -ffree-line-length-none -I$(VegasDir) -I$(ModuleDir) -J$(ModuleDir) $(LHAPDFflags) -Dcompiler=2 
endif


ifeq ($(Comp),ifort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O2 -vec-report0 -diag-disable remark
 else
    OptOpts  = -O0 -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv  -diag-disable remark 
 endif
 ifeq ($(UseLHAPDF),Yes)
    LHAPDFflags = -L$(MyLHADir) -lLHAPDF -DuseLHAPDF=1
 else
    LHAPDFflags = -DuseLHAPDF=0
 endif 
 fcomp = ifort -fpp $(OptOpts) -module $(ModuleDir) -I$(VegasDir) -Dcompiler=1 $(LHAPDFflags)
endif


# never use gcc with other than O0, optimization is buggy
ccomp = gcc -O0


# executable file
Exec = ./JHUGen


# main dependencies
MainDep = mod_Kinematics.F90 \
          mod_CrossSection.F90 \
          mod_CrossSection_TTBH.F90 \
          mod_CrossSection_BBBH.F90 \
          main.F90

# main objects
MainObj = $(ObjectDir)/mod_Kinematics.o \
          $(ObjectDir)/mod_CrossSection.o \
          $(ObjectDir)/mod_CrossSection_TTBH.o \
          $(ObjectDir)/mod_CrossSection_BBBH.o \
          $(ObjectDir)/main.o


# amplitude dependencies
AmpDep = mod_Misc.F90 \
          mod_Parameters.F90 \
          mod_TopDecay.F90 \
          mod_Graviton.F90 \
          mod_Higgs.F90 \
          mod_Zprime.F90 \
          mod_HiggsJJ.F90 \
          mod_HiggsJ.F90 \
          mod_VHiggs.F90 \
          mod_TTBHiggs.F90


# amplitude objects
AmpObj =  $(ObjectDir)/mod_Misc.o \
          $(ObjectDir)/mod_Parameters.o \
          $(ObjectDir)/mod_TopDecay.o \
          $(ObjectDir)/mod_Graviton.o \
          $(ObjectDir)/mod_Higgs.o \
          $(ObjectDir)/mod_Zprime.o \
          $(ObjectDir)/mod_HiggsJJ.o \
          $(ObjectDir)/mod_HiggsJ.o \
          $(ObjectDir)/mod_VHiggs.o \
          $(ObjectDir)/mod_TTBHiggs.o

VegasDep = $(VegasDir)/vegas.f
VegasObj = $(ObjectDir)/vegas.o

PDFDep   = $(PDFDir)/Cteq61Pdf.f \
           $(PDFDir)/mstwpdf.f \
           $(PDFDir)/NNPDFDriver.f

PDFObj   = $(ObjectDir)/Cteq61Pdf.o \
           $(ObjectDir)/mstwpdf.o \
           $(ObjectDir)/NNPDFDriver.o

PSGenDep = genps.c \
           boost.c
PSGenObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o


# MadGraphDep = h_emepemep.f \
#               h_emepemep_noint.f \
#               h_tamtaptamtap.f \
#               h_tamtaptamtap_noint.f \
#               switchmom.f
MadGraphDep = gg_ttbh.f \
              uub_ttbh.f \
              switchmom.f

MadGraphObj = $(ObjectDir)/gg_ttbh.o \
              $(ObjectDir)/uub_ttbh.o \
              $(ObjectDir)/switchmom.o \
              $(HOME)/lib/HELAS-3.0/libdhelas3.ifc90.a



JHUGenerator: $(PSGenObj) $(VegasObj) $(PDFObj) $(AmpObj) $(MainObj)
	@echo " linking"
ifeq ($(UseLHAPDF),Yes)
	echo " interfaced with LHAPDF"
else
	echo " using internal PDF sets" 
endif
	@echo " executable file is "$(Exec) "compiled with "$(Comp)
	@echo " "
	$(fcomp) -o $(Exec) $(MainObj) $(AmpObj) $(PSGenObj) $(VegasObj) $(PDFObj)



$(AmpObj): $(AmpDep) mod_Misc.F90 mod_Parameters.F90 $(Makefile)
	@echo " compiling amplitude routines"
	$(fcomp) -c $(AmpDep)
	mv *.o $(ObjectDir)


$(MainObj): $(MainDep) $(VegasDep) $(PDFDep) $(PSGenDep) mod_Parameters.F90 $(Makefile)
	@echo " compiling main routines"
	$(fcomp) -c $(MainDep)
	mv *.o $(ObjectDir)

$(VegasObj): $(VegasDep) $(Makefile)
	@echo " compiling vegas"
	$(fcomp) -c $(VegasDep)
	mv *.o $(ObjectDir)

$(PDFObj): $(PDFDep) $(Makefile)
	@echo " compiling pdfs"
	$(fcomp) -c $(PDFDep)
	mv *.o $(ObjectDir)

$(PSGenObj): $(PSGenDep)
	@echo " compiling phase space generator"
	$(ccomp) -c $(PSGenDep)
	mv *.o $(ObjectDir)


# $(MadGraphObj): $(MadGraphDep)
# 	@echo " compiling MG files"
# 	$(fcomp) -O0 -c $(MadGraphDep) -I/home/schulze/lib/HELAS-3.0/
# 	mv *.o $(ObjectDir)



clean:
	@echo " deleting executable, object and module files"
	rm -f objects/*.o
	rm -f modules/*.mod
	rm -f modules/*.f90
	rm -f *.mod
	rm -f $(Exec)

# supresses command calls
.SILENT:
