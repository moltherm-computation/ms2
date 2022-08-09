# CRAY.mk
# for CRAY XE6,XC40,... systems
# using Cray compiler
# prerequisite: module load PrgEnv-cray
# (should be named CRAYcray.mk)
#
F90              := ftn
F90ld            := $(F90)
F90mpi           := $(F90)
#OMPFLAGS         := 	# OpenMP is processed as default - no option required
OMPFLAGS          = -h omp
CPPFLAGS          = -DARCH=2 -DFORTRAN=2003
ifeq ($(OMP),1)
F90FLAGS_RELEASE  = -O3 -e m -rm -hvector3 -hpl=$(shell pwd)/pl
F90FLAGS_DEBUG    = -O0 -e m -g -hbounds
F90FLAGS_PROF     = -O3 -e m -g -rm
else
# noomp flag to disable OpenMP
F90FLAGS_RELEASE  = -O3 -e m -h noomp -rm -hvector3 -hpl=$(shell pwd)/pl
F90FLAGS_DEBUG    = -O0 -e m -h noomp -g -hbounds
F90FLAGS_PROF     = -O3 -e m -h noomp -g -rm
endif
PROG              = ms2_CRAY

LDFLAGS_RELEASE   = -O3 -rm -hwp -hpl=$(shell pwd)/pl
LDFLAGS_DEBUG     = -g 
LDFLAGS_PROF      = -O3 -g -rm

# remarks:
#  for F90FLAGS_RELEASE,LDFLAGS_RELEASE:
#   -hwp -hpl=/fullpath/builddir/PL.1  for whole program mode
