Here = $(PWD)
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
PDFDir = $(Here)/pdfs
VegasDir = $(Here)/vegas

# compiler options
Opt   = Yes
ifeq ($(Opt),Yes)
   IfortOpts  = -O2 -fpp -vec-report0
else
   IfortOpts  = -O0 -fpp -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv
endif
fcomp = gfort $(IfortOpts) -module $(ModuleDir) -I$(VegasDir)

# never use gcc with other than O0, optimization is buggy
ccomp = gcc -O0


# executable file
Exec = ./graviton


# main dependencies
MainDep = mod_Misc.f90 \
          mod_Parameters.f90 \
          mod_Kinematics.f90 \
          mod_Graviton.f90 \
          mod_Higgs.f90 \
          mod_Zprime.f90 \
          mod_CrossSection.f90 \
          main.f90

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
	@echo " executabel file is " $(Exec)
	@echo " "
	$(fcomp) -o $(Exec) $(MainObj) $(PSGenObj) $(VegasObj) $(PDFObj)


$(MainObj): $(MainDep) $(VegasDep) $(PDFDep)  $(PSGenDep) makefile
	@echo " compiling main routines"
	$(fcomp) -c $(MainDep)
	mv *.o $(ObjectDir)

$(VegasObj): $(VegasDep)
	@echo " compiling vegas"
	$(fcomp) -c $(VegasDep)
	mv *.o $(ObjectDir)

$(PDFObj): $(PDFDep)
	@echo " compiling pdfs"
	$(fcomp) -c $(PDFDep)
	mv *.o $(ObjectDir)

$(PSGenObj): $(PSGenDep)
	@echo " compiling phase space generator"
	$(ccomp) -c $(PSGenDep)
	mv *.o $(ObjectDir)



# supresses command calls
.SILENT:
