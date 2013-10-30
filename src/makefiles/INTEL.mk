F90              = ifort
OMPFLAGS         = -openmp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90  
F90FLAGS_RELEASE = -O3 -fpp -r8 -vec_report0 
F90FLAGS_DEBUG   = -O0 -fpp -r8 -vec_report0 -g 
F90FLAGS_PROF    = -O3 -fpp -r8 -vec_report0 -g 
PROG             = ms2_INTEL

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -g

