# --------------------- Macro-Defs -------------------------
F90 = gfortran
# --------------------- End-Macro-Defs ----------

# Link step
Forward:DIFF_SOLVER.o ODE.o Normal.o Random.o Forward.o
	$(F90) -o Forward DIFF_SOLVER.o ODE.o Normal.o Random.o Forward.o

# Compiling steps

DIFF_SOLVER.o:./DIFF_SOLVER.f95
	$(F90) -O3 -c ./DIFF_SOLVER.f95

ODE.o:./ODE.f95 DIFF_SOLVER.o
	$(F90) -O3 -c ./ODE.f95

Normal.o:./Normal.f95
	$(F90) -O3 -c ./Normal.f95

Random.o:./Random.f95 Normal.o
	$(F90) -O3 -c ./Random.f95

Forward.o:./Forward.f95 ODE.o Random.o
	$(F90) -O3 -c ./Forward.f95
