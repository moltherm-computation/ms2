F90              = ftn
F90ld            = ftn
F90mpi           = ftn
OMPFLAGS         = -fopenmp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90 
F90FLAGS_RELEASE = -O3 -ffree-line-length-none  -mtune=generic
F90FLAGS_DEBUG   = -O0 -ffree-line-length-none  -mtune=generic -g
F90FLAGS_PROF    = -O3 -ffree-line-length-none  -mtune=generic -pg
PROG             = ms2_CRAYgnu

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -pg

