Comp = gfort

ifeq ($(Comp),ifort)
fcomp = ifort -fpp -O2 -vec-report0 -Dcompiler=1
# fcomp = ifort -fpp -O0 -vec-report0 -Dcompiler=1
endif
ifeq ($(Comp),gfort)
fcomp = f95 -O3 -ffree-line-length-none -Dcompiler=2
# fcomp = f95 -O0 -ffree-line-length-none -Dcompiler=2 -fcheck=all

endif

ccomp = gcc


MCFM_Dep = WBFZZ_MCFM/src/WBFZZ/anomamp.f \
	   WBFZZ_MCFM/src/WBFZZ/qq_ZZqq.f \
           WBFZZ_MCFM/src/WBFZZ/setupzprops.f \
           WBFZZ_MCFM/src/WBFZZ/spinorcurr.f \
           WBFZZ_MCFM/src/WBFZZ/jzero.f \
           WBFZZ_MCFM/src/WBFZZ/jone.f \
           WBFZZ_MCFM/src/WBFZZ/jtwo.f \
           WBFZZ_MCFM/src/WBFZZ/jtwo3456.f \
           WBFZZ_MCFM/src/WBFZZ/ZZSingleres.f \
           WBFZZ_MCFM/src/WBFZZ/ZZHZZamp.f \
           WBFZZ_MCFM/src/WBFZZ/WWZZ.f \
           WBFZZ_MCFM/src/Need/cdotpr.f

# MCFM_Obj = WBFZZ_MCFM/src/WBFZZ/*.o \
#            WBFZZ_MCFM/src/Need/cdotpr.o
MCFM_Obj = anomamp.o cdotpr.o jtwo3456.o jtwo.o jzero.o qq_ZZqq.o setupzprops.o spinorcurr.o WWZZ.o ZZHZZamp.o ZZSingleres.o jone.o



Testprogram: mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_HiggsJ_MatEl.o mod_HiggsJJ_MatEl.o mod_VHiggs_MatEl.o mod_TTBH_MatEl.o mod_TH_MatEl.o vegas.o NNPDFDriver.o testprogram.F90 \
	checkWBF_SM.dat checkWBF_1-8.dat checkSBF_SM.dat checkSBF_1-4.dat checkHJ_SM.dat checkZH_SM.dat MCFMforVBF
	@echo " "
	@echo " compiling and linking testprogram.F90 with "$(Comp)
	$(fcomp) -o testF testprogram.F90 -lm vegas.o NNPDFDriver.o mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_HiggsJ_MatEl.o mod_HiggsJJ_MatEl.o mod_VHiggs_MatEl.o mod_TTBH_MatEl.o mod_TH_MatEl.o $(MCFM_Obj)
	@echo " "
	@echo " compiling and linking testprogram.c with gcc"
	$(ccomp) -o testC testprogram.c NNPDFDriver.o vegas.o mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_TTBH_MatEl.o mod_TH_MatEl.o -lm -lgfortran
	@echo " "


NNPDFDriver.o: ./pdfs/NNPDFDriver.f
	@echo " "
	@echo " compiling NNPDFDriver.f with "$(Comp)
	$(fcomp) -c ./pdfs/NNPDFDriver.f
	
	
vegas.o: ./vegas.f
	@echo " "
	@echo " compiling vegas.f with "$(Comp)
	$(fcomp) -c ./vegas.f


mod_Higgs_MatEl.o: mod_Higgs_MatEl.F90 includeVars.F90 variables.F90
	@echo " "
	@echo " compiling mod_Higgs_MatEl.F90 includeVars.F90 with "$(Comp)
	$(fcomp) -c mod_Higgs_MatEl.F90  


mod_Zprime_MatEl.o: mod_Zprime_MatEl.F90 includeVars.F90 variables.F90
	@echo " "
	@echo " compiling mod_Zprime_MatEl.F90 includeVars.F90 with "$(Comp)
	$(fcomp) -c mod_Zprime_MatEl.F90  


mod_Graviton_MatEl.o: mod_Graviton_MatEl.F90 includeVars.F90 variables.F90
	@echo " "
	@echo " compiling mod_Graviton_MatEl.F90 includeVars.F90 with "$(Comp)
	$(fcomp) -c mod_Graviton_MatEl.F90  

mod_HiggsJ_MatEl.o: mod_HiggsJ_MatEl.F90 variables.F90
	@echo " "
	@echo " compiling mod_HiggsJ_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJ_MatEl.F90

mod_HiggsJJ_MatEl.o: mod_HiggsJJ_MatEl.F90 variables.F90
	@echo " "
	@echo " compiling mod_HiggsJJ_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJJ_MatEl.F90

mod_VHiggs_MatEl.o: mod_VHiggs_MatEl.F90 variables.F90
	@echo " "
	@echo " compiling mod_VHiggs_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_VHiggs_MatEl.F90 

mod_TTBH_MatEl.o: mod_TTBH_MatEl.F90 variables.F90
	@echo " "
	@echo " compiling mod_TTBH_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_TTBH_MatEl.F90  

mod_TH_MatEl.o: mod_TH_MatEl.F90 variables.F90
	@echo " "
	@echo " compiling mod_TH_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_TH_MatEl.F90  

MCFMforVBF: $(MCFM_Dep)
	@echo " "
	@echo " compiling MCFM WBF files with "$(Comp)
	$(fcomp) -c -I./WBFZZ_MCFM/src/Inc/ $(MCFM_Dep)
# 
# 	@echo " compiling MCFM WBF files with ifort!"
# 	ifort  -O0 -implicitnone -zero -check bounds -check pointer -warn interfaces -ftrapuv  -diag-disable remark -debug extended -g -traceback -fpe0 -check uninit  -c -I./WBFZZ_MCFM/src/Inc/ $(MCFM_Dep)
# 	ifort  -O2  -c -I./WBFZZ_MCFM/src/Inc/ $(MCFM_Dep)

clean:
	@echo " deleting object files"
	rm -f *.o *.mod
	rm -f testF testC



# supresses command calls
.SILENT:
