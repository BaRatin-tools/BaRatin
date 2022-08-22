#========================================================================
# Makefile to build BaM for LINUX based on MiniDMSL
#========================================================================
#

# Define the Fortran Compiler and its options
FC  = gfortran
FLAGS = -fcheck=all -Wall -O3 # DEBUG: -fcheck=all -Wall -g -O0

# Define directories
DMSL_DIR = ../../miniDMSL/src/
BMSL_DIR = ../../BMSL/src/
BaRatin_DIR = ../src/

# Define DMSL files
DMSL = \
	$(DMSL_DIR)kinds_dmsl_kit.f90              \
	$(DMSL_DIR)types_dmsl_kit.f90              \
	$(DMSL_DIR)utilities_dmsl_kit.f90          \
	$(DMSL_DIR)mt19937-64.f95                  \
	$(DMSL_DIR)uniran1_minidmsl_mod.f90        \
	$(DMSL_DIR)linalg_dmsl_kit.f90             \
	$(DMSL_DIR)numerix_dmsl_kit.f90

# Define BMSL files
BMSL = \
	$(BMSL_DIR)Dates_tools.f90                   \
	$(BMSL_DIR)Aggregation_tools.f90             \
	$(BMSL_DIR)EmpiricalStats_tools.f90          \
	$(BMSL_DIR)TimeSeries_tools.f90              \
	$(BMSL_DIR)Geodesy_tools.f90                 \
	$(BMSL_DIR)DataRW_tools.f90                  \
	$(BMSL_DIR)Distribution_tools.f90            \
	$(BMSL_DIR)MCMC_tools.f90                    \
	$(BMSL_DIR)MCMCStrategy_tools.f90            \
	$(BMSL_DIR)BayesianEstimation_tools.f90 


# Define BaRatin files
BaRatin = \
	$(BaRatin_DIR)RatingCurve_tools.f90        \
	$(BaRatin_DIR)BaRatin_tools.f90            \
	$(BaRatin_DIR)Main.f90 

all:build_exe

build_exe: compile_baratin
	$(FC)  $(FLAGS) $(DMSL) $(BMSL) $(BaRatin) -o BaRatin

compile_baratin: compile_bmsl
	$(FC) $(FLAGS) -c $(BaRatin)

compile_bmsl: compile_dmsl
	$(FC) $(FLAGS) -c $(BMSL)

compile_dmsl:
	$(FC) $(FLAGS) -c $(DMSL)





