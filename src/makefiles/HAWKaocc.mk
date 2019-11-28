# for CRAY XE6,XC40,... systems
# using Portland Group Inc. (PGI) compiler
# prerequisite: module load PrgEnv-pgi
#
F90              = mpif90
F90ld            = mpif90
F90mpi           = mpif90
OMPFLAGS         = -mp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90 -D_PGF
F90FLAGS_RELEASE = -r8 -Bstatic -finline-aggressive -flto -march=znver2
F90FLAGS_DEBUG   = -g -Mdwarf2 -r8 -Minfo
#F90FLAGS_PROF    = -Mpfi -Mpfo
PROG             = ms2_HAWKaocc

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
#LDFLAGS_PROF     = 

