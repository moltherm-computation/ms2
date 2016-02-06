F90              = ftn
F90ld            = ftn
F90mpi           = ftn
OMPFLAGS         = -openmp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90  
F90FLAGS_RELEASE = -O3 -fpp -r8 -qopt-report=5
F90FLAGS_DEBUG   = -O0 -fpp -r8 -qopt-report=0 -g 
F90FLAGS_PROF    = -O3 -fpp -r8 -qopt-report=0 -g 
PROG             = ms2_CRAYintel

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -g

