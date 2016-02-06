F90              = ftn
F90ld            = ftn
F90mpi           = ftn
OMPFLAGS         = -h omp
CPPFLAGS         = -DARCH=2 -DFORTRAN=2003 
F90FLAGS_RELEASE = -O3 -e m -h noomp -rm -hvector3
F90FLAGS_DEBUG   = -O0 -e m -h noomp -g
F90FLAGS_PROF    = -O3 -e m -h noomp -g -rm
PROG             = ms2_CRAY

LDFLAGS_RELEASE  = -O3 -rm
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -g -rm

