F90              = gfortran
F90ld            = gfortran
F90mpi           = mpif90
OMPFLAGS         = -fopenmp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90
# GCC 10.X:                                                   -fallow-argument-mismatch
F90FLAGS_RELEASE = -O3 -ffree-line-length-none -mtune=generic 
F90FLAGS_DEBUG   = -O0 -ffree-line-length-none -mtune=generic -g -fcheck=all
F90FLAGS_PROF    = -O3 -ffree-line-length-none -mtune=generic -pg
PROG             = ms2_GNU

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -pg

