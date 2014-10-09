F90              = ftn
OMPFLAGS         = -h omp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90 
F90FLAGS_RELEASE = -O3 -e m -h noomp -rm -hvector3
F90FLAGS_DEBUG   = -O0 -e m -h noomp -g
F90FLAGS_PROF    = -O3 -e m -h noomp -g -rm
PROG             = ms2_CRAY

LDFLAGS_RELEASE  = -O3 -rm
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -g -rm

