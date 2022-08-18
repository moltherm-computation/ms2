F90              = gfortran
F90ld            = gfortran
F90mpi           = mpif90

# Initialize variable
VERSIONDEPENDENDOPTION =

ifeq ($(MPI), 1) # if MPI is active

    # Get complete version, e.g. 5.4.0
    GCCVERSION = $(shell $(F90mpi) -dumpfullversion)

    # Get the first (major) version number, e.g. 5 from 5.4.0
    MAJORGCCVERSION = $(shell echo $(GCCVERSION) | cut -f1 -d.)

    # version is 10.* or 11.*?
    ifneq ($(filter $(MAJORGCCVERSION),10 11),)

        $(info Version of $(F90mpi) is $(MAJORGCCVERSION).*: add flag '-fallow-argument-mismatch' and hide all warnings)
        VERSIONDEPENDENDOPTION += -fallow-argument-mismatch -w
        # replace errors with warnings and hide ALL warnings
    endif
endif

OMPFLAGS         = -fopenmp
CPPFLAGS         = -DARCH=2 -DFORTRAN=90

F90FLAGS_RELEASE = -O3 -ffree-line-length-none -mtune=generic

# add version dependent flags if not empty, 'strip' removes redundant white space
ifneq ($(VERSIONDEPENDENDOPTION),)
     F90FLAGS_RELEASE += $(strip $(VERSIONDEPENDENDOPTION))
endif

F90FLAGS_DEBUG   = -O0 -ffree-line-length-none -mtune=generic -g -fcheck=all
F90FLAGS_PROF    = -O3 -ffree-line-length-none -mtune=generic -pg
PROG             = ms2_GNU

LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -g 
LDFLAGS_PROF     = -O3 -pg

