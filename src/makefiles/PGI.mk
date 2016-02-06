F90              = pgf90
F90ld            = pgf90
F90mpi           = mpif90
OMPFLAGS         = -mp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90 -D_PGF 
F90FLAGS_RELEASE = -C -O3 
F90FLAGS_DEBUG   = -C -O0 -g 
F90FLAGS_PROF    = -C -O3 -g 
PROG             = ms2_PGI

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -g

