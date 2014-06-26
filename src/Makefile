########################################
# Makefile for ms2
########################################
COMPILER = GNU
CONFIG = makefiles/$(COMPILER).mk
include $(CONFIG)
TARGET    = RELEASE
MPI       = 0
OMP       = 0
SP        = 0
PROG      = ms2
TRANS     = 0
ifeq ($(MPI),1)
  ifneq ($(COMPILER),CRAY)
     F90 = mpif90
  endif
  CPPFLAGS += -DMPI_VER=1
  PROG     :=  $(PROG)_$(COMPILER)_mpi
else
  CPPFLAGS += -DMPI_VER=0
  PROG     :=  $(PROG)_$(COMPILER)
endif

F90FLAGS = $(F90FLAGS_$(TARGET))
LDFLAGS  = $(LDFLAGS_$(TARGET))

ifeq ($(OMP),1)
CPPFLAGS += -DENABLE_OMP
F90FLAGS += $(OMPFLAGS)
LDFLAGS  += $(OMPFLAGS)
PROG     := $(PROG)_omp
endif

ifeq ($(SP),1)
CPPFLAGS += -DSINGLEPRECISION
LDFLAGS  += -DSINGLEPRECISION
PROG     := $(PROG)_SP
endif

ifeq ($(TRANS),1)
CPPFLAGS += -DTRANS
PROG     := $(PROG)_trans
endif

SRCS =  ms2.F90 \
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
	@make clean
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(OBJS) $(MODS)

.SUFFIXES: .F90 .mod .o

%.mod %.o: %.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<

########################################
# Dependencies for modules
#
ms2.o: ms2_simulation.mod ms2_global.mod

ms2_simulation.mod: ms2_ensemble.mod ms2_accumulator.mod ms2_global.mod \
  ms2_stopwatch.mod

ms2_ensemble.mod: ms2_component.mod ms2_interaction.mod \
  ms2_accumulator.mod ms2_global.mod

ms2_component.mod: ms2_molecule.mod ms2_site.mod \
  ms2_accumulator.mod ms2_global.mod

ms2_molecule.mod: ms2_site.mod ms2_global.mod

ms2_site.mod: ms2_global.mod

ms2_interaction.mod: ms2_potential.mod ms2_component.mod ms2_global.mod

ms2_potential.mod: ms2_site.mod ms2_global.mod

ms2_stopwatch.mod: ms2_global.mod
