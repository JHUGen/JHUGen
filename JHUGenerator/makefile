Here = $(PWD)
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas

# compiler options, Comp=ifort/gfort
Opt  = Yes
Comp = ifort



ifeq ($(Comp),gfort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O3
 else
    OptOpts  = -O0
 endif
 fcomp = f95 $(OptOpts) -ffree-line-length-none -I$(VegasDir) -J$(ModuleDir) -Dcompiler=2
endif


ifeq ($(Comp),ifort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O2 -vec-report0
 else
    OptOpts  = -O0 -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv
 endif
 fcomp = ifort -fpp $(OptOpts) -module $(ModuleDir) -I$(VegasDir) -Dcompiler=1
endif


# never use gcc with other than O0, optimization is buggy
ccomp = gcc -O0


# executable file
Exec = ./JHUGen


# main dependencies
MainDep = mod_Misc.F90 \
          mod_Parameters.F90 \
          mod_Kinematics.F90 \
          mod_Graviton.F90 \
          mod_Higgs.F90 \
          mod_Zprime.F90 \
          mod_CrossSection.F90 \
          main.F90

# main objects
MainObj = $(ObjectDir)/mod_Misc.o \
          $(ObjectDir)/mod_Parameters.o \
          $(ObjectDir)/mod_Kinematics.o \
          $(ObjectDir)/mod_Graviton.o \
          $(ObjectDir)/mod_Higgs.o \
          $(ObjectDir)/mod_Zprime.o \
          $(ObjectDir)/mod_CrossSection.o \
          $(ObjectDir)/main.o


VegasDep = $(VegasDir)/vegas.f
VegasObj = $(ObjectDir)/vegas.o

PDFDep   = $(PDFDir)/Cteq61Pdf.f \
           $(PDFDir)/mstwpdf.f

PDFObj   = $(ObjectDir)/Cteq61Pdf.o \
           $(ObjectDir)/mstwpdf.o

PSGenDep = genps.c \
           boost.c
PSGenObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o



JHUGenerator: $(PSGenObj) $(VegasObj) $(PDFObj) $(MainObj)
	@echo " linking"
	@echo " executable file is "$(Exec) "compiled with "$(Comp)
	@echo " "
	$(fcomp) -o $(Exec) $(MainObj) $(PSGenObj) $(VegasObj) $(PDFObj)


$(MainObj): $(MainDep) $(VegasDep) $(PDFDep)  $(PSGenDep) makefile
	@echo " compiling main routines"
	$(fcomp) -c $(MainDep)
	mv *.o $(ObjectDir)

$(VegasObj): $(VegasDep) makefile
	@echo " compiling vegas"
	$(fcomp) -c $(VegasDep)
	mv *.o $(ObjectDir)

$(PDFObj): $(PDFDep) makefile
	@echo " compiling pdfs"
	$(fcomp) -c $(PDFDep)
	mv *.o $(ObjectDir)

$(PSGenObj): $(PSGenDep)
	@echo " compiling phase space generator"
	$(ccomp) -c $(PSGenDep)
	mv *.o $(ObjectDir)



clean:
	@echo " deleting executable, object and module files"
	rm -f objects/*.o
	rm -f modules/*.mod
	rm -f *.mod
	rm -f $(Exec)

# supresses command calls
.SILENT:
