F90              = sxf90
#F90              = sxf03
F90ld            = sxf90
#F90ld            = sxf03
F90mpi           = sxmpif90
#F90mpi           = sxmpif03
OMPFLAGS         = 
CPPFLAGS         = -DARCH=3 -DFORTRAN=90 -f5
F90FLAGS_RELEASE = -C hopt -pi auto -Wf'-pvctl vchg loopfusion loopchg outerstrip outerunroll'
##F90FLAGS_DEBUG   = -g -EP -ftrace
F90FLAGS_DEBUG   = -C debug -EP
F90FLAGS_PROF    = 
PROG             = ms2_SX

LDFLAGS_RELEASE  = 
LDFLAGS_DEBUG    = 
LDFLAGS_PROF     = 

