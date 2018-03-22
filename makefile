# --------------------- Macro-Defs -------------------------
F90 = gfortran
# --------------------- End-Macro-Defs ----------

# Link step
MetropolisTest:MVSTAT.o DIFF_SOLVER.o PsiModule.o Pdf.o GetModule.o ODE.o GradientModule.o Normal.o Random.o Reading.o Sampling.o Prior.o Metropolis.o
	$(F90) -o MetropolisTest MVSTAT.o DIFF_SOLVER.o PsiModule.o Pdf.o GetModule.o ODE.o GradientModule.o Normal.o Random.o Reading.o Sampling.o Prior.o Metropolis.o

# Compiling steps

MVSTAT.o:./MVSTAT.f90
	$(F90) -O3 -c ./MVSTAT.f90

DIFF_SOLVER.o:./DIFF_SOLVER.f95
	$(F90) -O3 -c ./DIFF_SOLVER.f95

PsiModule.o:./PsiModule.f90
	$(F90) -O3 -c ./PsiModule.f90

GetModule.o:./GetModule.f95
	$(F90) -O3 -c ./GetModule.f95

Normal.o:./Normal.f95
	$(F90) -O3 -c ./Normal.f95
Random.o:./Random.f95 Normal.o
	$(F90) -O3 -c ./Random.f95

Reading.o:./Reading.f95
	$(F90) -O3 -c ./Reading.f95

ODE.o:./ODE.f95 DIFF_SOLVER.o
	$(F90) -O3 -c ./ODE.f95

Pdf.o:./Pdf.f95 MVSTAT.o ODE.o
	$(F90) -O3 -c ./Pdf.f95


GradientModule.o:./GradientModule.f95 PsiModule.o GetModule.o Pdf.o ODE.o
	$(F90) -O3 -c ./GradientModule.f95

Sampling.o:./Sampling.f95 GradientModule.o ODE.o Normal.o Random.o
	$(F90) -O3 -c ./Sampling.f95

Prior.o:./Prior.f95 Sampling.o Normal.o ODE.o Random.o
	$(F90) -O3 -c ./Prior.f95

Metropolis.o:./Metropolis.f95 Sampling.o Prior.o Reading.o
	$(F90) -O3 -c ./Metropolis.f95

# Make clean enabler
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
