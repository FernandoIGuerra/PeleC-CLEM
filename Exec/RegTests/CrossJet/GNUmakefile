# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Reference case for CLEM validation: plain PeleC, no CLEM subgrid model.
CLEM_MODEL = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
PELE_CVODE_FORCE_YCORDER = FALSE
PELE_USE_MAGMA = FALSE
PELE_COMPILE_AJACOBIAN = TRUE
Eos_Model := Fuego
Chemistry_Model := BurkeDryer
Transport_Model := Simple

# GNU Make
Bpack := ./Make.package
Blocs := .
PELE_HOME := ../../..
include $(PELE_HOME)/Exec/Make.PeleC
