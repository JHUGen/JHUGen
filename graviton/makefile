Here = $(PWD)
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas

# compiler options
Opt  = Yes
Comp = gfort



ifeq ($(Comp),gfort)
 ifeq ($(Opt),Yes)
    OptOpts  = -O3
 else
    OptOpts  = -O0
 endif
# fcomp = f95 -cpp $(OptOpts) -ffree-line-length-none -J$(ModuleDir) -I$(VegasDir) -Dcompiler=2
 fcomp = f95 $(OptOpts) -ffree-line-length-none -I$(VegasDir) -Dcompiler=2
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
Exec = ./graviton


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

PDFDep   = $(PDFDir)/Cteq61Pdf.f
PDFObj   = $(ObjectDir)/Cteq61Pdf.o

PSGenDep = genps.c \
           boost.c
PSGenObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o



graviton: $(PSGenObj) $(VegasObj) $(PDFObj) $(MainObj)
	@echo " linking"
	@echo " executabel file is "$(Exec) "compiled with "$(Comp)
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
	@echo " deleting object and module files"
	rm -f objects/*
	rm -f modules/*
	rm -f *.mod

# supresses command calls
.SILENT:
