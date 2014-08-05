Comp = gfort

ifeq ($(Comp),ifort)
fcomp = ifort -fpp -O2 -vec-report0 -Dcompiler=1
# fcomp = ifort -fpp -O0 -vec-report0 -Dcompiler=1
endif
ifeq ($(Comp),gfort)
fcomp = f95 -O3 -ffree-line-length-none -Dcompiler=2
endif

ccomp = gcc



Testprogram: mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_HiggsJJ_MatEl.o testprogram.F90 \
	checkWBF_SM.dat checkWBF_1-8.dat checkSBF_SM.dat checkSBF_1-4.dat
	@echo " "
	@echo " compiling and linking testprogram.F90 with "$(Comp)
	$(fcomp) -o testF testprogram.F90 -lm mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o mod_HiggsJJ_MatEl.o
	@echo " "
	@echo " compiling and linking testprogram.c with gcc"
	$(ccomp) -o testC testprogram.c -lm -lgfortran mod_Higgs_MatEl.o mod_Zprime_MatEl.o mod_Graviton_MatEl.o
	@echo " "




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

mod_HiggsJJ_MatEl.o: mod_HiggsJJ_MatEl.F90 variables.F90
	@echo " "
	@echo " compiling mod_HiggsJJ_MatEl.F90 with "$(Comp)
	$(fcomp) -c mod_HiggsJJ_MatEl.F90  

clean:
	@echo " deleting object files"
	rm -f *.o *.mod
	rm -f testF testC



# supresses command calls
.SILENT:
