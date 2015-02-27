Comp = gfort

ifeq ($(Comp),ifort)
fcomp = ifort -fpp -O2 -vec-report0 -Dcompiler=1
# fcomp = ifort -fpp -O0 -vec-report0 -Dcompiler=1
endif
ifeq ($(Comp),gfort)
fcomp = f95 -O3 -ffree-line-length-none -Dcompiler=2
endif

ccomp = gcc



Testprogram: mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_HiggsJ_MatEl.o mod_HiggsJJ_MatEl.o mod_VHiggs_MatEl.o mod_TTBH_MatEl.o NNPDFDriver.o testprogram.F90 \
	checkWBF_SM.dat checkWBF_1-8.dat checkSBF_SM.dat checkSBF_1-4.dat checkHJ_SM.dat checkZH_SM.dat
	@echo " "
	@echo " compiling and linking testprogram.F90 with "$(Comp)
	$(fcomp) -o testF testprogram.F90 -lm NNPDFDriver.o mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_HiggsJ_MatEl.o mod_HiggsJJ_MatEl.o mod_VHiggs_MatEl.o mod_TTBH_MatEl.o
	@echo " "
	@echo " compiling and linking testprogram.c with gcc"
	$(ccomp) -o testC testprogram.c NNPDFDriver.o mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_TTBH_MatEl.o  -lm -lgfortran
	@echo " "



NNPDFDriver.o: ./pdfs/NNPDFDriver.f
	@echo " "
	@echo " compiling NNPDFDriver.f with "$(Comp)
	$(fcomp) -c ./pdfs/NNPDFDriver.f


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
		
clean:
	@echo " deleting object files"
	rm -f *.o *.mod
	rm -f testF testC



# supresses command calls
.SILENT:
