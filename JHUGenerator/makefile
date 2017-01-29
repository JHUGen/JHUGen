Here = $(PWD)
Makefile=makefile
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
MELADir = $(Here)/../JHUGenMELAv2/MELA
MELAFortranDir = $(MELADir)/fortran
MELADataDir = $(MELADir)/data
MELALibDir = $(MELADataDir)/slc6_amd64_gcc530
MCFMLib = mcfm_703 # libmcfm_703.so
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas

# compiler options, Comp='ifort' or 'gfort'
Opt  = Yes
Comp = gfort

# link pdfs via LHA library ('Yes' or 'No')
UseLHAPDF=No
# remember to export
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}

# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
# this one in terms of the environment variable should work in most cases
MyLHADir=${LHAPDF_DATA_PATH}/../../lib/

#link MCFM libraries from MELA
linkMELA = No


ifeq ($(UseLHAPDF),Yes)
   LHAPDFflags = -L$(MyLHADir) -lLHAPDF -DuseLHAPDF=1
else
   LHAPDFflags = -DuseLHAPDF=0
endif

MELAflags = -DlinkMELA=0
ifeq ($(linkMELA),Yes)
 ifeq ($(Comp),gfort)
    MELAflags = -DlinkMELA=1 -L$(MELALibDir) -l$(MCFMLib)
 else
    linkMELA = No
 endif
endif

ifeq ($(Comp),gfort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O3
 else
    OptOpts  = -O0 -fcheck=all
 endif
 fcomp = gfortran $(OptOpts) -ffree-line-length-none -I$(VegasDir) -I$(ModuleDir) -J$(ModuleDir) $(LHAPDFflags) $(MELAflags) -Dcompiler=2
endif


ifeq ($(Comp),ifort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O2 -vec-report0 -diag-disable remark
 else
    OptOpts  = -O0 -implicitnone -check bounds -check pointer -warn interfaces -ftrapuv  -diag-disable remark -debug extended -g -traceback -fpe0 -check uninit
 endif
 fcomp = ifort -fpp $(OptOpts) -module $(ModuleDir) -I$(VegasDir) -Dcompiler=1 $(LHAPDFflags) $(MELAflags)
endif


# never use gcc with other than O0, optimization is buggy
ccomp = gcc -O0


# executable file
Exec = ./JHUGen


# kinematics dependencies
KinDep = mod_Phasespace.F90 \
         mod_Kinematics.F90

# kinematics objects
KinObj = $(ObjectDir)/mod_Phasespace.o \
         $(ObjectDir)/mod_Kinematics.o

ifeq ($(linkMELA),Yes)
    KinDep += mod_MCFMWrapper.F90
    KinObj += $(ObjectDir)/mod_MCFMWrapper.o
endif

# main dependencies
MainDep = mod_CrossSection.F90 \
          mod_CrossSection_HJJ.F90 \
          mod_CrossSection_TTBH.F90 \
          mod_CrossSection_BBBH.F90 \
          mod_CrossSection_TH.F90 \
          mod_PMZZ.F90 \
          main.F90

# main objects
MainObj = $(ObjectDir)/mod_CrossSection.o \
          $(ObjectDir)/mod_CrossSection_HJJ.o \
          $(ObjectDir)/mod_CrossSection_TTBH.o \
          $(ObjectDir)/mod_CrossSection_BBBH.o \
          $(ObjectDir)/mod_CrossSection_TH.o \
          $(ObjectDir)/mod_PMZZ.o \
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
          mod_TTBHiggs.F90 \
          mod_THiggs.F90


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
          $(ObjectDir)/mod_TTBHiggs.o \
          $(ObjectDir)/mod_THiggs.o

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

CPSDep = ./CPS/CALLING_cpHTO.f

CPSObj = $(ObjectDir)/CALLING_cpHTO.o


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



JHUGenerator: $(CPSObj) $(PSGenObj) $(VegasObj) $(PDFObj) $(AmpObj) $(KinObj) $(MainObj)
	@echo " linking"
ifeq ($(UseLHAPDF),Yes)
	echo " -> interfaced with LHAPDF"
else
	echo " -> using internal PDF sets"
endif
	@echo " executable file "$(Exec)" is compiled with "$(Comp)
	@echo " "
	$(fcomp) -o $(Exec) $(MainObj) $(KinObj) $(AmpObj) $(PSGenObj) $(VegasObj) $(PDFObj) $(CPSObj) $(LHAPDFflags) $(MELAflags)



$(AmpObj): $(AmpDep) mod_Misc.F90 mod_Parameters.F90 $(Makefile)
	@echo " compiling amplitude routines"
	$(fcomp) -c $(AmpDep)
	mv *.o $(ObjectDir)

$(KinObj): $(KinDep) $(VegasDep) $(PDFDep) $(PSGenDep) mod_Parameters.F90 $(Makefile)
	@echo " compiling kinematics routines"
ifeq ($(linkMELA),Yes)
	echo " -> interfaced with JHUGenMELA and MCFM"
	ln -sf $(MELADataDir)/input.DAT ./
	ln -sf $(MELADataDir)/process.DAT ./
	ln -sf $(MELADataDir)/hto_output.dat ./
	ln -sf $(MELADataDir)/ffwarn.dat ./
	ln -sf $(MELADataDir)/ffperm5.dat ./
	ln -sf $(MELADataDir)/fferr.dat ./
	ln -sf $(MELADataDir)/dm_parameters.DAT ./
	ln -sf $(MELADataDir)/br.sm1 ./
	ln -sf $(MELADataDir)/br.sm2 ./
endif
	$(fcomp) -c $(KinDep)
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

$(CPSObj): $(CPSDep) $(Makefile)
	@echo " compiling complex-pole scheme for the spin-0 propagator"
	$(fcomp) -c $(CPSDep)
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
ifeq ($(linkMELA),Yes)
	rm -f input.DAT
	rm -f process.DAT
	rm -f hto_output.dat
	rm -f ffwarn.dat
	rm -f ffperm5.dat
	rm -f fferr.dat
	rm -f dm_parameters.DAT
	rm -f br.sm1
	rm -f br.sm2
endif
	rm -f objects/*.o
	rm -f modules/*.mod
	rm -f modules/*.f90
	rm -f *.mod
	rm -f *.o
	rm -f $(Exec)

# supresses command calls
.SILENT:
