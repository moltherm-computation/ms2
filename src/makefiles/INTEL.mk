F90              = ifort
F90ld            = ifort
F90mpi           = mpiifort
OMPFLAGS         = -qopenmp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90
F90FLAGS_RELEASE = -O3 -fpp -r8 -qopt-report=5 
F90FLAGS_DEBUG   = -O0 -fpp -r8 -qopt-report=0 -g -debug all -check all
F90FLAGS_PROF    = -O3 -fpp -r8 -qopt-report=0 -g 
PROG             = ms2_INTEL

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -g

# remarks:
#  for F90FLAGS_RELEASE,LDFLAGS_RELEASE:
#   -fast (= -ipo -O3 -no-prec-div -static -fp-model fast=2 -xHost)
#         might raise problems due to implied -static
#         with -ipo the compile time might be a little faster, but link and overall time will raise
#  
#  for AMD platforms default settings or -xHost might not give best optimizations...

