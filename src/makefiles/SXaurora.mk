F90              = nfort
F90ld            = nfort
F90mpi           = mpinfort
OMPFLAGS         =
CPPFLAGS         = -DARCH=3 -DFORTRAN=90
#                  -O3: segmentation faults
F90FLAGS_RELEASE = -O2 -report-all
#F90FLAGS_DEBUG   = -g -traceback #-ftrace
F90FLAGS_DEBUG   = -O0 -report-all -traceback -fcheck=bounds -Werror
F90FLAGS_PROF    =
PROG             = ms2_AURORA

LDFLAGS_RELEASE  =
LDFLAGS_DEBUG    = -traceback
LDFLAGS_PROF     =
