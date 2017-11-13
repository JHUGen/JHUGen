Here = $(PWD)
Makefile=makefile
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
MELADir = $(Here)/../JHUGenMELAv2/MELA
#MELADir = $(Here)/../../Ulascan/JHUGen/JHUGenMELAv2/MELA
MELAFortranDir = $(MELADir)/fortran
MELADataDir = $(MELADir)/data
MELALibDir = $(MELADataDir)/slc6_amd64_gcc530
MCFMLib = mcfm_705 # libmcfm_705.so
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas
CollierDir = $(Here)/COLLIER-1.2/

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
# remember to export
#          LD_LIBRARY_PATH=/.../MELA/data/$SCRAM_ARCH/:${LD_LIBRARY_PATH}
linkMELA = No


ifeq ($(UseLHAPDF),Yes)
   LHAPDFflags = -L$(MyLHADir) -lLHAPDF -I$(MyLHADir)/../include -DuseLHAPDF=1 -lstdc++
else
   LHAPDFflags = -DuseLHAPDF=0
endif


#include Collier Loop Integral Library
linkCollierLib = Yes

ifeq ($(linkCollierLib),Yes)
#   CollierLibflags = -DuseCollier=1 -L$(CollierDir) -I$(CollierDir)/modules -lcollier
   CollierLibflags = -DuseCollier=1 -L$(CollierDir) -I$(CollierDir) -lcollier
   #COLLIERflags = -L$(COLLIERDir) -lcollier
else
   CollierLibflags = -DuseCollier=0
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
 fcomp = gfortran $(OptOpts) -ffree-line-length-none -I$(VegasDir) -I$(ModuleDir) -J$(ModuleDir) $(LHAPDFflags) $(MELAflags) $(CollierLibflags) -Dcompiler=2
endif


ifeq ($(Comp),ifort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O2 -vec-report0 -diag-disable remark
 else
    OptOpts  = -O0 -implicitnone -check bounds -check pointer -warn interfaces -ftrapuv  -diag-disable remark -debug extended -g -traceback -fpe0 -check uninit
 endif
 fcomp = ifort -fpp $(OptOpts) -module $(ModuleDir) -I$(VegasDir) -Dcompiler=1 $(LHAPDFflags) $(MELAflags) $(CollierLibflags)
endif


# never use gcc with other than O0, optimization is buggy
ccomp = gcc -O0


# executable file
Exec = ./JHUGen


# kinematics dependencies
KinDep = mod_HashCollection.F90 \
         mod_Phasespace.F90 \
         mod_Kinematics.F90 \
         mod_Kinematics_VH.F90

# kinematics objects
KinObj = $(ObjectDir)/mod_HashCollection.o \
         $(ObjectDir)/mod_Phasespace.o \
         $(ObjectDir)/mod_Kinematics.o \
         $(ObjectDir)/mod_Kinematics_VH.o

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
          mod_CrossSection_VH.F90 \
          mod_PMZZ.F90 \
          main.F90

# main objects
MainObj = $(ObjectDir)/mod_CrossSection.o \
          $(ObjectDir)/mod_CrossSection_HJJ.o \
          $(ObjectDir)/mod_CrossSection_TTBH.o \
          $(ObjectDir)/mod_CrossSection_BBBH.o \
          $(ObjectDir)/mod_CrossSection_TH.o \
          $(ObjectDir)/mod_CrossSection_VH.o \
          $(ObjectDir)/mod_PMZZ.o \
          $(ObjectDir)/main.o

# auxiliary dependencies
AuxDep = mod_Misc.F90 \
         mod_Parameters.F90 \
         mod_VH_aux.F90 \

AuxObj =  $(ObjectDir)/mod_Misc.o \
          $(ObjectDir)/mod_Parameters.o \
          $(ObjectDir)/mod_VH_aux.o \

ggboxVHAmpDep = \
         gg_VH_Boxes/mod_ffbH_Box7_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box8_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box8_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box8_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box8_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box10_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box10_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box10_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box10_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box12_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box12_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box12_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box12_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box7_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box8_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box8_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box9_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box9_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box10_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box10_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box11_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box11_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box12_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box12_mpp.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box8_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box8_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box10_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box10_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box12_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box12_pmp.F90

ggboxVHAmpObj = \
          $(ObjectDir)/mod_ffbH_Box7_mmm.o \
          $(ObjectDir)/mod_ffbH_Box7_ppm.o \
          $(ObjectDir)/mod_ffbH_Box7_mmp.o \
          $(ObjectDir)/mod_ffbH_Box7_ppp.o \
          $(ObjectDir)/mod_ffbH_Box8_mmm.o \
          $(ObjectDir)/mod_ffbH_Box8_ppm.o \
          $(ObjectDir)/mod_ffbH_Box8_mmp.o \
          $(ObjectDir)/mod_ffbH_Box8_ppp.o \
          $(ObjectDir)/mod_ffbH_Box9_mmm.o \
          $(ObjectDir)/mod_ffbH_Box9_ppm.o \
          $(ObjectDir)/mod_ffbH_Box9_mmp.o \
          $(ObjectDir)/mod_ffbH_Box9_ppp.o \
          $(ObjectDir)/mod_ffbH_Box10_mmm.o \
          $(ObjectDir)/mod_ffbH_Box10_ppm.o \
          $(ObjectDir)/mod_ffbH_Box10_mmp.o \
          $(ObjectDir)/mod_ffbH_Box10_ppp.o \
          $(ObjectDir)/mod_ffbH_Box11_mmm.o \
          $(ObjectDir)/mod_ffbH_Box11_ppm.o \
          $(ObjectDir)/mod_ffbH_Box11_mmp.o \
          $(ObjectDir)/mod_ffbH_Box11_ppp.o \
          $(ObjectDir)/mod_ffbH_Box12_mmm.o \
          $(ObjectDir)/mod_ffbH_Box12_ppm.o \
          $(ObjectDir)/mod_ffbH_Box12_mmp.o \
          $(ObjectDir)/mod_ffbH_Box12_ppp.o \
          $(ObjectDir)/mod_ffbH_Box7_mpm.o \
          $(ObjectDir)/mod_ffbH_Box7_mpp.o \
          $(ObjectDir)/mod_ffbH_Box8_mpm.o \
          $(ObjectDir)/mod_ffbH_Box8_mpp.o \
          $(ObjectDir)/mod_ffbH_Box9_mpm.o \
          $(ObjectDir)/mod_ffbH_Box9_mpp.o \
          $(ObjectDir)/mod_ffbH_Box10_mpm.o \
          $(ObjectDir)/mod_ffbH_Box10_mpp.o \
          $(ObjectDir)/mod_ffbH_Box11_mpm.o \
          $(ObjectDir)/mod_ffbH_Box11_mpp.o \
          $(ObjectDir)/mod_ffbH_Box12_mpm.o \
          $(ObjectDir)/mod_ffbH_Box12_mpp.o \
          $(ObjectDir)/mod_ffbH_Box7_pmm.o \
          $(ObjectDir)/mod_ffbH_Box7_pmp.o \
          $(ObjectDir)/mod_ffbH_Box8_pmm.o \
          $(ObjectDir)/mod_ffbH_Box8_pmp.o \
          $(ObjectDir)/mod_ffbH_Box9_pmm.o \
          $(ObjectDir)/mod_ffbH_Box9_pmp.o \
          $(ObjectDir)/mod_ffbH_Box10_pmm.o \
          $(ObjectDir)/mod_ffbH_Box10_pmp.o \
          $(ObjectDir)/mod_ffbH_Box11_pmm.o \
          $(ObjectDir)/mod_ffbH_Box11_pmp.o \
          $(ObjectDir)/mod_ffbH_Box12_pmm.o \
          $(ObjectDir)/mod_ffbH_Box12_pmp.o

# amplitude dependencies
AmpDep =  mod_Parameters.F90 \
          mod_Misc.F90 \
          mod_TopDecay.F90 \
          mod_Graviton.F90 \
          mod_Higgs.F90 \
          mod_Zprime.F90 \
          mod_HiggsJJ.F90 \
          mod_HiggsJ.F90 \
          mod_VHiggs.F90 \
          mod_VH_LO.F90 \
          mod_VH_gg.F90 \
          mod_VH_real.F90 \
          mod_VH_qg.F90 \
          mod_VH_dipole.F90 \
          mod_TTBHiggs.F90 \
          mod_THiggs.F90


# amplitude objects
AmpObj =  $(ObjectDir)/mod_Parameters.o \
          $(ObjectDir)/mod_Misc.o \
          $(ObjectDir)/mod_TopDecay.o \
          $(ObjectDir)/mod_Graviton.o \
          $(ObjectDir)/mod_Higgs.o \
          $(ObjectDir)/mod_Zprime.o \
          $(ObjectDir)/mod_HiggsJJ.o \
          $(ObjectDir)/mod_HiggsJ.o \
          $(ObjectDir)/mod_VHiggs.o \
          $(ObjectDir)/mod_VH_LO.o \
          $(ObjectDir)/mod_VH_gg.o \
          $(ObjectDir)/mod_VH_real.o \
          $(ObjectDir)/mod_VH_qg.o \
          $(ObjectDir)/mod_VH_dipole.o \
          $(ObjectDir)/mod_TTBHiggs.o \
          $(ObjectDir)/mod_THiggs.o

ggVHAmpObj = $(ObjectDir)/ampBox_ggVH.o
ggVHAmpDep = ampBox_ggVH.F90
          

VegasDep = $(VegasDir)/vegas.f
VegasObj = $(ObjectDir)/vegas.o

PDFDep   = $(PDFDir)/Cteq61Pdf.f \
           $(PDFDir)/mstwpdf.f \
           $(PDFDir)/NNPDFDriver.f

PDFObj   = $(ObjectDir)/Cteq61Pdf.o \
           $(ObjectDir)/mstwpdf.o \
           $(ObjectDir)/NNPDFDriver.o

PSGenC = genps.c \
         boost.c

PSGenDep = $(PSGenC) \
           genps.h \
           boost.h

PSGenObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o

CPSDep = ./CPS/CALLING_cpHTO.f

CPSObj = $(ObjectDir)/CALLING_cpHTO.o

LHAPDFInterfaceC = LHAPDFversion.cc

LHAPDFInterfaceDep = $(LHAPDFInterfaceC) \
                     LHAPDFversion.h

ifeq ($(UseLHAPDF),Yes)
LHAPDFInterfaceObj = $(ObjectDir)/LHAPDFversion.o
else
LHAPDFInterfaceObj =
endif


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




JHUGenerator: colliermodules $(CPSObj) $(PSGenObj) $(VegasObj) $(PDFObj) $(AuxObj) $(ggboxVHAmpObj) $(AmpObj) $(ggVHAmpObj) $(KinObj) $(MainObj)
	@echo " linking"
ifeq ($(UseLHAPDF),Yes)
	echo " -> interfaced with LHAPDF"
else
	echo " -> using internal PDF sets"
endif
	@echo " executable file "$(Exec)" is compiled with "$(Comp)
	@echo " "
	$(fcomp) -o $(Exec) $(MainObj) $(AuxObj) $(ggboxVHAmpObj) $(KinObj) $(AmpObj) $(ggVHAmpObj) $(PSGenObj) $(VegasObj) $(PDFObj) $(CPSObj) $(LHAPDFflags) $(MELAflags) 

$(AuxObj): $(AuxDep) $(Makefile)
	@echo " compiling misc. and aux. tools"
	$(fcomp) -c $(AuxDep)
	mv *.o $(ObjectDir)

$(ggboxVHAmpObj): $(ggboxVHAmpDep) $(Makefile)
	@echo " compiling ggVH box amplitude"
	$(fcomp) -c $(ggboxVHAmpDep)
	mv *.o $(ObjectDir)

$(AmpObj): $(AmpDep) mod_Misc.F90 mod_Parameters.F90 $(Makefile)
	@echo " compiling amplitude routines"
	$(fcomp) -c $(AmpDep)
	mv *.o $(ObjectDir)

$(ggVHAmpObj): $(ggVHAmpDep) $(Makefile)
	@echo " compiling ggVH box amplitude"
	$(fcomp) -c $(ggVHAmpDep)
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

$(PSGenObj): $(PSGenDep) $(Makefile)
	@echo " compiling phase space generator"
	$(ccomp) -c $(PSGenC)
	mv *.o $(ObjectDir)

ifeq ($(UseLHAPDF),Yes)
$(LHAPDFInterfaceObj) : $(LHAPDFInterfaceDep) $(Makefile)
	@echo " compiling LHAPDF interface"
	$(ccomp) $(LHAPDFflags) -c $(LHAPDFInterfaceC)
	mv *.o $(ObjectDir)
endif

# $(MadGraphObj): $(MadGraphDep)
# 	@echo " compiling MG files"
# 	$(fcomp) -O0 -c $(MadGraphDep) -I/home/schulze/lib/HELAS-3.0/
# 	mv *.o $(ObjectDir)

colliermodules:
#	@echo " compiling COLLIER-1.2"
#	FC=$(CompFullName) overrides FC specified in that makefile
#	(cd $(CollierDir); $(MAKE) FC=$(Comp))
	@echo " moving COLLIER .mod and .dylib to modules directory"
	(cp $(CollierDir)/*.mod $(ModuleDir)/)
	(cp $(CollierDir)/*.dylib $(ModuleDir)/)
#	@echo " moving Collier .o to objects directory"
#	(cp $(CollierDir)/*.o $(objects)/)

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
