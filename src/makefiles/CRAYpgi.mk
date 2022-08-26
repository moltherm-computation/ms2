# for CRAY XE6,XC40,... systems
# using Portland Group Inc. (PGI) compiler
# prerequisite: module load PrgEnv-pgi
#
F90              = ftn
F90ld            = ftn
F90mpi           = ftn
OMPFLAGS         = -mp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90 -D_PGF -V
F90FLAGS_RELEASE = -fastsse -r8 -Bstatic -Mipa=fast,inline
F90FLAGS_DEBUG   = -g -Mdwarf2 -r8 -Minfo
#F90FLAGS_PROF    = -Mpfi -Mpfo
PROG             = ms2_CRAYpgi

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
#LDFLAGS_PROF     = 

