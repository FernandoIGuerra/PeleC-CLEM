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

# CLEM subgrid model. TRUE builds the CLEM-active executable used by
# inputs-case1-nonreacting-CLEM.inp; CLEM_NLEM is the compile-time number of
# LEM elements per LES cell (ClemIndex.H, default 12).
#
# The plain-PeleC reference for the same case is the SAME executable name, so
# the two builds overwrite each other. Build the reference with
#   make CLEM_MODEL=FALSE -j20
# and keep the binary aside, or use a separate copy of this directory.
CLEM_MODEL = TRUE
CLEM_NLEM  = 12

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
