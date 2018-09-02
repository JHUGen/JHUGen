Here = $(PWD)
Makefile=makefile
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
MELADir = $(Here)/../JHUGenMELA/MELA
#MELADir = $(Here)/../../Ulascan/JHUGen/JHUGenMELA/MELA
MELAFortranDir = $(MELADir)/fortran
MELADataDir = $(MELADir)/data
MELALibDir = $(MELADataDir)/slc6_amd64_gcc530
MCFMLib = mcfm_705 # libmcfm_705.so
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas
CollierDir = $(Here)/../lib/COLLIER-1.1/

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
linkMELA = Yes

ifeq ($(UseLHAPDF),Yes)
   LHAPDFflags = -L$(MyLHADir) -lLHAPDF -I$(MyLHADir)/../include -DuseLHAPDF=1 -lstdc++
else
   LHAPDFflags = -DuseLHAPDF=0
endif


#include Collier Loop Integral Library
linkCollierLib = No

ifeq ($(linkMELA),Yes)
   linkCollierLib = Yes
   CollierLibflags = -DuseCollier=1 -I$(MELADir)/COLLIER -lcollier
else
   ifeq ($(linkCollierLib),Yes)
      CollierLibflags = -DuseCollier=1 -L$(CollierDir) -I$(CollierDir)/modules -lcollier
   else
      CollierLibflags = -DuseCollier=0
   endif
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
 CompFullName = gfortran
 ifeq ($(Opt),Yes)
    OptOpts  = -O3
 else
    OptOpts  = -O0 -fcheck=all
 endif
 fcomp = gfortran $(OptOpts) -ffree-line-length-none -I$(VegasDir) -I$(ModuleDir) -J$(ModuleDir) $(LHAPDFflags) $(MELAflags) $(CollierLibflags) -Dcompiler=2
endif


ifeq ($(Comp),ifort)
 CompFullName = ifort
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


# main dependencies
MainDep = mod_CrossSection.F90 \
          mod_CrossSection_HJJ.F90 \
          mod_CrossSection_TTBH.F90 \
          mod_CrossSection_BBBH.F90 \
          mod_CrossSection_TH.F90 \
          mod_CrossSection_VH.F90 \
          mod_PMZZ.F90 \
          mod_CommandLine.F90 \
          main.F90

ifeq ($(linkCollierLib),Yes)
		 MainDep += mod_CrossSection_HH.F90
endif

# main objects
MainObj = $(ObjectDir)/mod_CrossSection.o \
          $(ObjectDir)/mod_CrossSection_HJJ.o \
          $(ObjectDir)/mod_CrossSection_TTBH.o \
          $(ObjectDir)/mod_CrossSection_BBBH.o \
          $(ObjectDir)/mod_CrossSection_TH.o \
          $(ObjectDir)/mod_CrossSection_VH.o \
          $(ObjectDir)/mod_PMZZ.o \
          $(ObjectDir)/mod_CommandLine.o \
          $(ObjectDir)/main.o

ifeq ($(linkCollierLib),Yes)
		 MainObj += $(ObjectDir)/mod_CrossSection_HH.o
endif

AuxDep = mod_Parameters.F90 \
         mod_Misc.F90 \
         mod_VH_aux.F90 \
         mod_Phasespace.F90 \
         mod_HashCollection.F90

KimDep = mod_Kinematics.F90 \
         mod_Kinematics_VH.F90

ifeq ($(linkCollierLib),Yes)
ggVHAmpDep = \
         gg_VH_Boxes/mod_ffbH_Box7_8_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_8_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_8_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_8_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_10_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_10_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_10_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_10_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_12_mmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_12_ppm.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_12_mmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_12_ppp.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_8_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box7_8_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box9_10_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box9_10_mpp.F90 \
		 gg_VH_Boxes/mod_ffbH_Box11_12_mpm.F90 \
		 gg_VH_Boxes/mod_ffbH_Box11_12_mpp.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_8_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box7_8_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_10_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box9_10_pmp.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_12_pmm.F90 \
         gg_VH_Boxes/mod_ffbH_Box11_12_pmp.F90

ggVH5AmpDep = \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_mmm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_ppm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_mmp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_ppp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box9_10_mmm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box9_10_ppm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box9_10_mmp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box9_10_ppp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box11_12_mmm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box11_12_ppm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box11_12_mmp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box11_12_ppp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_mpm.F90 \
		 gg_VH5_Boxes/mod_ffbH5_Box7_8_mpp.F90 \
		 gg_VH5_Boxes/mod_ffbH5_Box9_10_mpm.F90 \
		 gg_VH5_Boxes/mod_ffbH5_Box9_10_mpp.F90 \
		 gg_VH5_Boxes/mod_ffbH5_Box11_12_mpm.F90 \
		 gg_VH5_Boxes/mod_ffbH5_Box11_12_mpp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_pmm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box7_8_pmp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box9_10_pmm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box9_10_pmp.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box11_12_pmm.F90 \
         gg_VH5_Boxes/mod_ffbH5_Box11_12_pmp.F90

ggHHAmpDep = \
		 gg_HH_Boxes/mod_gg_box_HH_1_mm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_1_mp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_1_pm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_1_pp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_2_mm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_2_mp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_2_pm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_2_pp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_3_mm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_3_mp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_3_pm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_3_pp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_4_mm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_4_mp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_4_pm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_4_pp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_5_mm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_5_mp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_5_pm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_5_pp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_6_mm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_6_mp.F90 \
         gg_HH_Boxes/mod_gg_box_HH_6_pm.F90 \
         gg_HH_Boxes/mod_gg_box_HH_6_pp.F90

ggHH5AmpDep = \
         gg_HH5_Boxes/mod_gg_box_HH5_1_mm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_1_mp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_1_pm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_1_pp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_2_mm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_2_mp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_2_pm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_2_pp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_3_mm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_3_mp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_3_pm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_3_pp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_4_mm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_4_mp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_4_pm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_4_pp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_5_mm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_5_mp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_5_pm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_5_pp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_6_mm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_6_mp.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_6_pm.F90 \
         gg_HH5_Boxes/mod_gg_box_HH5_6_pp.F90

ggVHAmpObj = \
          $(ObjectDir)/mod_ffbH_Box7_8_mmm.o \
          $(ObjectDir)/mod_ffbH_Box7_8_ppm.o \
          $(ObjectDir)/mod_ffbH_Box7_8_mmp.o \
          $(ObjectDir)/mod_ffbH_Box7_8_ppp.o \
          $(ObjectDir)/mod_ffbH_Box9_10_mmm.o \
          $(ObjectDir)/mod_ffbH_Box9_10_ppm.o \
          $(ObjectDir)/mod_ffbH_Box9_10_mmp.o \
          $(ObjectDir)/mod_ffbH_Box9_10_ppp.o \
          $(ObjectDir)/mod_ffbH_Box11_12_mmm.o \
          $(ObjectDir)/mod_ffbH_Box11_12_ppm.o \
          $(ObjectDir)/mod_ffbH_Box11_12_mmp.o \
          $(ObjectDir)/mod_ffbH_Box11_12_ppp.o \
          $(ObjectDir)/mod_ffbH_Box7_8_mpm.o \
          $(ObjectDir)/mod_ffbH_Box7_8_mpp.o \
          $(ObjectDir)/mod_ffbH_Box9_10_mpm.o \
          $(ObjectDir)/mod_ffbH_Box9_10_mpp.o \
          $(ObjectDir)/mod_ffbH_Box11_12_mpm.o \
          $(ObjectDir)/mod_ffbH_Box11_12_mpp.o \
          $(ObjectDir)/mod_ffbH_Box7_8_pmm.o \
          $(ObjectDir)/mod_ffbH_Box7_8_pmp.o \
          $(ObjectDir)/mod_ffbH_Box9_10_pmm.o \
          $(ObjectDir)/mod_ffbH_Box9_10_pmp.o \
          $(ObjectDir)/mod_ffbH_Box11_12_pmm.o \
          $(ObjectDir)/mod_ffbH_Box11_12_pmp.o

ggVH5AmpObj = \
          $(ObjectDir)/mod_ffbH5_Box7_8_mmm.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_ppm.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_mmp.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_ppp.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_mmm.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_ppm.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_mmp.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_ppp.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_mmm.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_ppm.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_mmp.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_ppp.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_mpm.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_mpp.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_mpm.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_mpp.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_mpm.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_mpp.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_pmm.o \
          $(ObjectDir)/mod_ffbH5_Box7_8_pmp.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_pmm.o \
          $(ObjectDir)/mod_ffbH5_Box9_10_pmp.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_pmm.o \
          $(ObjectDir)/mod_ffbH5_Box11_12_pmp.o

ggHHAmpObj = \
          $(ObjectDir)/mod_gg_box_HH_1_mm.o \
          $(ObjectDir)/mod_gg_box_HH_1_mp.o \
          $(ObjectDir)/mod_gg_box_HH_1_pm.o \
          $(ObjectDir)/mod_gg_box_HH_1_pp.o \
          $(ObjectDir)/mod_gg_box_HH_2_mm.o \
          $(ObjectDir)/mod_gg_box_HH_2_mp.o \
          $(ObjectDir)/mod_gg_box_HH_2_pm.o \
          $(ObjectDir)/mod_gg_box_HH_2_pp.o \
          $(ObjectDir)/mod_gg_box_HH_3_mm.o \
          $(ObjectDir)/mod_gg_box_HH_3_mp.o \
          $(ObjectDir)/mod_gg_box_HH_3_pm.o \
          $(ObjectDir)/mod_gg_box_HH_3_pp.o \
          $(ObjectDir)/mod_gg_box_HH_4_mm.o \
          $(ObjectDir)/mod_gg_box_HH_4_mp.o \
          $(ObjectDir)/mod_gg_box_HH_4_pm.o \
          $(ObjectDir)/mod_gg_box_HH_4_pp.o \
          $(ObjectDir)/mod_gg_box_HH_5_mm.o \
          $(ObjectDir)/mod_gg_box_HH_5_mp.o \
          $(ObjectDir)/mod_gg_box_HH_5_pm.o \
          $(ObjectDir)/mod_gg_box_HH_5_pp.o \
          $(ObjectDir)/mod_gg_box_HH_6_mm.o \
          $(ObjectDir)/mod_gg_box_HH_6_mp.o \
          $(ObjectDir)/mod_gg_box_HH_6_pm.o \
          $(ObjectDir)/mod_gg_box_HH_6_pp.o

ggHH5AmpObj = \
          $(ObjectDir)/mod_gg_box_HH5_1_mm.o \
          $(ObjectDir)/mod_gg_box_HH5_1_mp.o \
          $(ObjectDir)/mod_gg_box_HH5_1_pm.o \
          $(ObjectDir)/mod_gg_box_HH5_1_pp.o \
          $(ObjectDir)/mod_gg_box_HH5_2_mm.o \
          $(ObjectDir)/mod_gg_box_HH5_2_mp.o \
          $(ObjectDir)/mod_gg_box_HH5_2_pm.o \
          $(ObjectDir)/mod_gg_box_HH5_2_pp.o \
          $(ObjectDir)/mod_gg_box_HH5_3_mm.o \
          $(ObjectDir)/mod_gg_box_HH5_3_mp.o \
          $(ObjectDir)/mod_gg_box_HH5_3_pm.o \
          $(ObjectDir)/mod_gg_box_HH5_3_pp.o \
          $(ObjectDir)/mod_gg_box_HH5_4_mm.o \
          $(ObjectDir)/mod_gg_box_HH5_4_mp.o \
          $(ObjectDir)/mod_gg_box_HH5_4_pm.o \
          $(ObjectDir)/mod_gg_box_HH5_4_pp.o \
          $(ObjectDir)/mod_gg_box_HH5_5_mm.o \
          $(ObjectDir)/mod_gg_box_HH5_5_mp.o \
          $(ObjectDir)/mod_gg_box_HH5_5_pm.o \
          $(ObjectDir)/mod_gg_box_HH5_5_pp.o \
          $(ObjectDir)/mod_gg_box_HH5_6_mm.o \
          $(ObjectDir)/mod_gg_box_HH5_6_mp.o \
          $(ObjectDir)/mod_gg_box_HH5_6_pm.o \
          $(ObjectDir)/mod_gg_box_HH5_6_pp.o
else
ggVHAmpDep = 
ggVH5AmpDep = 
ggHHAmpDep = 
ggHH5AmpDep = 
ggVHAmpObj = 
ggVH5AmpObj = 
ggHHAmpObj = 
ggHH5AmpObj = 
endif

# amplitude dependencies
AmpDep = mod_TopDecay.F90 \
         mod_Graviton.F90 \
         mod_Higgs.F90 \
         mod_Zprime.F90 \
         mod_VHiggs.F90 \
         mod_HiggsJJ.F90 \
         mod_HiggsJ.F90 \
         mod_VH_LO.F90 \
         mod_VH_real.F90 \
         mod_VH_qg.F90 \
         mod_VH_dipole.F90 \
         mod_TTBHiggs.F90 \
         mod_THiggs.F90

AuxObj =  $(ObjectDir)/mod_HashCollection.o \
          $(ObjectDir)/mod_Parameters.o \
          $(ObjectDir)/mod_Misc.o \
          $(ObjectDir)/mod_VH_aux.o \
          $(ObjectDir)/mod_Phasespace.o

KimObj = $(ObjectDir)/mod_Kinematics.o \
          $(ObjectDir)/mod_Kinematics_VH.o

ifeq ($(linkMELA),Yes)
    AuxDep += mod_MCFMWrapper.F90
    AuxObj += $(ObjectDir)/mod_MCFMWrapper.o
endif


AmpObj =  $(ObjectDir)/mod_TopDecay.o \
          $(ObjectDir)/mod_Graviton.o \
          $(ObjectDir)/mod_Higgs.o \
          $(ObjectDir)/mod_Zprime.o \
          $(ObjectDir)/mod_VHiggs.o \
          $(ObjectDir)/mod_HiggsJJ.o \
          $(ObjectDir)/mod_HiggsJ.o \
          $(ObjectDir)/mod_VH_LO.o \
          $(ObjectDir)/mod_VH_real.o \
          $(ObjectDir)/mod_VH_qg.o \
          $(ObjectDir)/mod_VH_dipole.o \
          $(ObjectDir)/mod_TTBHiggs.o \
          $(ObjectDir)/mod_THiggs.o

ifeq ($(linkCollierLib),Yes)
		 KimDep += mod_Kinematics_HH.F90
		 AmpDep += mod_VH_virtual.F90
         AmpDep += mod_VH_gg.F90
         AmpDep += mod_HH.F90
         AmpObj += $(ObjectDir)/mod_VH_gg.o
         AmpObj += $(ObjectDir)/mod_VH_virtual.o
         AmpObj += $(ObjectDir)/mod_HH.o
		 KimObj += $(ObjectDir)/mod_Kinematics_HH.o
endif    

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



JHUGenerator: $(CPSObj) $(PSGenObj) $(VegasObj) $(PDFObj) $(AuxObj) $(ggVHAmpObj) $(LHAPDFInterfaceObj) $(ggVH5AmpObj) $(ggHHAmpObj) $(ggHH5AmpObj) $(AmpObj) $(KimObj) $(MainObj)
	@echo " linkMELA = "$(linkMELA)
	@echo " linkCollierLib = "$(linkCollierLib)
	@echo " linking"
ifeq ($(UseLHAPDF),Yes)
	echo " -> interfaced with LHAPDF"
else
	echo " -> using internal PDF sets"
endif
	@echo " executable file "$(Exec)" is compiled with "$(Comp)
	@echo " "
	$(fcomp) -o $(Exec) $(MainObj) $(PSGenObj) $(AuxObj) $(AmpObj) $(KimObj) $(ggVHAmpObj) $(LHAPDFInterfaceObj) $(ggVH5AmpObj) $(ggHHAmpObj) $(ggHH5AmpObj) $(VegasObj) $(PDFObj) $(CPSObj) $(MCFM_Obj) $(LHAPDFflags)

$(AuxObj): $(AuxDep) $(Makefile)
	@echo " compiling misc. and aux. tools"
	$(fcomp) -c $(AuxDep)
	mv *.o $(ObjectDir)
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
	$(fcomp) -c $(AuxDep)
	mv *.o $(ObjectDir)

$(ggVHAmpObj): $(ggVHAmpDep) $(AuxObj) $(Makefile)
	@echo " compiling ggVH box amplitudes"
	$(fcomp) -c $(ggVHAmpDep)
	mv *.o $(ObjectDir)

$(ggVH5AmpObj): $(ggVH5AmpDep) $(AuxObj) $(Makefile)
	@echo " compiling ggVH5 box amplitudes"
	$(fcomp) -c $(ggVH5AmpDep)
	mv *.o $(ObjectDir)

$(ggHHAmpObj): $(ggHHAmpDep) $(AuxObj) $(Makefile)
	@echo " compiling ggHH box amplitudes"
	$(fcomp) -c $(ggHHAmpDep)
	mv *.o $(ObjectDir)

$(ggHH5AmpObj): $(ggHH5AmpDep) $(AuxObj) $(Makefile)
	@echo " compiling ggHH5 box amplitudes"
	$(fcomp) -c $(ggHH5AmpDep)
	mv *.o $(ObjectDir)

$(KimObj): $(KimDep) $(Makefile)
	@echo " compiling Kinematics dependencies"
	$(fcomp) -c $(KimDep)
	mv *.o $(ObjectDir)

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
	rm -f objects/*.a
	rm -f ff/*.o
	rm -f ff/*.a
	rm -f ql/*.o
	rm -f ql/*.a
	rm -f modules/*.mod
	rm -f modules/*.f90
	rm -f *.mod
	rm -f *.o
	rm -f $(Exec)

# supresses command calls
.SILENT:
