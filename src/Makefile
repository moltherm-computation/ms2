########################################
# Makefile for ms2
########################################
MPI=0
FVM=1

ifneq ($(MPI),0)
  ifeq ($(FVM),0)
    F90 = mpif90
    CPPFLAGS = -DARCH=2 -DFORTRAN=90 -DMPI_VER=1 -DFVM_VER=0 -DPAR_PROF
    #removed '-D_PGF'
    #symbol 'PAR_PROF' added by Hendrik Adorf (ITWM), 17.03.2009
    PROG =  ms2mpi
  else
    F90 = pgf90 
    CPPFLAGS = -DARCH=2 -DFORTRAN=90 -DMPI_VER=0 -DFVM_VER=0 -D_PGF
    PROG =  ms2
  endif
endif

ifeq ($(MPI),0)
  ifneq ($(FVM),0)
    F90 = ifort
    # use at least intel/10.1/64
    CPPFLAGS = -DARCH=2 -DFORTRAN=2003 -DMPI_VER=0 -DFVM_VER=1 #-DPAR_PROF -DPAR_DEBUG
    #removed '-D_PGF'
    FVMOBJS = /u/herc/adorf/lib/FVM-Fortran/libfvmf2003.o \
              /u/herc/adorf/lib/FVM-Fortran/libfvmf2003.a \
              /u/herc/adorf/lib/FVM-Fortran/fvmf2003extensions.o \
              /u/herc/adorf/lib/FVM-Fortran/fvmf2003extensions.a
    PROG =  ms2fvmf03
  else
    F90 = pgf90
    CPPFLAGS = -DARCH=2 -DFORTRAN=90 -DMPI_VER=0 -DFVM_VER=0 -D_PGF
    PROG =  ms2
  endif
endif

F90FLAGS = -g
#removed '-C'; Hendrik Adorf (ITWM), 17.03.2009
LDFLAGS =

SRCS =  ms2.F90 \
        ms2_profiler.F90 \
        ms2_stopwatch.F90 \
        ms2_global.F90 \
        ms2_simulation.F90 \
        ms2_ensemble.F90 \
        ms2_component.F90 \
        ms2_molecule.F90 \
        ms2_site.F90 \
        ms2_interaction.F90 \
        ms2_potential.F90 \
        ms2_accumulator.F90

OBJS =  $(foreach SRC,$(SRCS),$(basename $(SRC)).o)

MODS =  $(foreach SRC,$(SRCS),$(basename $(SRC)).mod)

########################################
# Rules to build the target
#
default: all

all: $(PROG)

$(PROG): $(OBJS)
ifneq ($(FVM),0)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(OBJS) $(FVMOBJS) \
	-L/u/herc/adorf/lib/PV-4D -libverbs12 -lPv4dVM4 \
	-lPv4dDaemon -lcrypto -lstdc++ \
	-Wl,-rpath=/u/herc/adorf/lib:/opt/cluster/Intel/cce/9.1.047/lib
else
	$(F90) $(LDFLAGS) -o $@ $(OBJS)
endif

clean:
	rm -f $(OBJS) $(MODS)

.SUFFIXES: .F90 .mod .o

%.mod %.o: %.F90
ifneq ($(FVM),0)
	$(F90) $(CPPFLAGS) $(F90FLAGS) \
	-I/u/herc/adorf/include/FVM-Fortran -c $<
else
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<
endif

########################################
# Dependencies for modules
#
ms2.o: ms2_simulation.mod ms2_global.mod

ms2_simulation.mod: ms2_ensemble.mod ms2_accumulator.mod ms2_global.mod \
  ms2_stopwatch.mod ms2_profiler.mod

ms2_ensemble.mod: ms2_component.mod ms2_interaction.mod \
  ms2_accumulator.mod ms2_global.mod ms2_profiler.mod

ms2_component.mod: ms2_molecule.mod ms2_site.mod \
  ms2_accumulator.mod ms2_global.mod ms2_profiler.mod

ms2_molecule.mod: ms2_site.mod ms2_global.mod ms2_profiler.mod

ms2_site.mod: ms2_global.mod

ms2_interaction.mod: ms2_potential.mod ms2_component.mod ms2_global.mod

ms2_potential.mod: ms2_site.mod ms2_global.mod

ms2_stopwatch.mod: ms2_global.mod

ms2_global.mod: ms2_profiler.mod
