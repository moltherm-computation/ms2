########################################
# Makefile for ms2
########################################
#
#---------------------------------------
#           see makefile/<COMPILER>.mk
COMPILER ?= GNU
#
TARGET   ?= RELEASE
#TARGET   ?= DEBUG
TRANS    ?= 0
OSMOP    ?= 0
HBOND    ?= 0
SPME     ?= 0
MPI      ?= 0
OMP      ?= 0
SP       ?= 0
#---------------------------------------
#
PROG     = ms2
PROG     :=  $(PROG)_$(COMPILER)
#
# including architecture specific definitions
CONFIG := makefiles/$(COMPILER).mk
include $(CONFIG)
#
# print some info
$(info COMPILER: $(COMPILER))
$(info TARGET: $(TARGET))
ifeq ($(TRANS),1)
$(info TRANSPORT enabled)
endif
ifeq ($(OSMOP),1)
$(info osmoticPressure enabled)
endif
ifeq ($(OSMOP),2)
$(info osmoticPressure with profiling enabled)
endif
ifeq ($(HBOND),1)
$(info H-Bonding statistics enabled)
endif
ifeq ($(SPME),1)
$(info Smooth Particle Mesh Ewald enabled)
endif
ifeq ($(MPI),1)
$(info MPI enabled)
endif
ifeq ($(OMP),1)
$(info OMP enabled)
endif
ifeq ($(SP),1)
$(info compiling Single Precision)
endif
#

F90FLAGS = $(F90FLAGS_$(TARGET))
LDFLAGS  = $(LDFLAGS_$(TARGET))

ifeq ($(TARGET),DEBUG)
  PROG     := $(PROG)_dbg
endif

ifeq ($(MPI),1)
  F90      = $(F90mpi)
  F90ld    = $(F90mpi)
  CPPFLAGS += -DMPI_VER=1
  PROG     :=  $(PROG)_mpi
else
  CPPFLAGS += -DMPI_VER=0
endif

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

ifeq ($(OSMOP),1)
  CPPFLAGS += -DOSMOP=1
  PROG     := $(PROG)_osmop1
endif

ifeq ($(OSMOP),2)
  CPPFLAGS += -DOSMOP=2
  PROG     := $(PROG)_osmop2
endif

ifeq ($(HBOND),1)
  CPPFLAGS += -DHBOND=1
  PROG     := $(PROG)_hbond
endif

ifeq ($(SPME),1)
  CPPFLAGS += -DSPME=1
  PROG     := $(PROG)_spme
endif

SRCS =  ms2.F90 \
        ms2_stopwatch.F90 \
        ms2_global.F90 \
        ms2_simulation.F90 \
        ms2_ensemble.F90 \
        ms2_component.F90 \
        math_types.F90 \
        ms2_molecule.F90 \
        ms2_unit.F90 \
        ms2_site.F90 \
        ms2_idf.F90 \
        ms2_interaction.F90 \
        ms2_potential.F90 \
        ms2_accumulator.F90

OBJS =  $(foreach SRC,$(SRCS),$(basename $(SRC)).o)

MODS =  $(foreach SRC,$(SRCS),$(basename $(SRC)).mod)

########################################
# Rules to build the target
#
.PHONY: clean

default: all
	@make clean

all: $(PROG)

ifeq ($(SPME),1)
$(PROG): $(OBJS)
	$(F90ld) $(LDFLAGS) -o $@ $(OBJS) -lfftw3 -lm 
endif
ifneq ($(SPME),1)
$(PROG): $(OBJS)
	$(F90ld) $(LDFLAGS) -o $@ $(OBJS)
endif

clean:
	rm -f $(OBJS) $(MODS)
	#rm -f *.optrpt

.SUFFIXES: .F90 .mod .o
ifeq ($(COMPILER),$(filter $(COMPILER), INTEL CRAYintel))
ifeq ($(TARGET),RELEASE)
# reduce optimization level from standard for Intel O3 to O2
# to prevent error 'segmentation violation signal raised'
ms2_component.mod ms2_component.o: ms2_molecule.mod ms2_site.mod \
                                   ms2_accumulator.mod ms2_global.mod
	$(F90) $(CPPFLAGS) $(F90FLAGS) -O2 -c ms2_component.F90
endif
endif

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
  ms2_accumulator.mod ms2_global.mod math_types.mod

math_types.mod: ms2_global.mod

ms2_molecule.mod: ms2_site.mod ms2_unit.mod ms2_idf.mod ms2_global.mod

ms2_unit.mod: ms2_site.mod ms2_global.mod

ms2_site.mod: ms2_global.mod

ms2_idf.mod:  ms2_global.mod

ms2_interaction.mod: ms2_potential.mod ms2_component.mod ms2_global.mod

ms2_potential.mod: ms2_site.mod ms2_idf.mod ms2_global.mod ms2_molecule.mod  math_types.mod

ms2_accumulator.mod: ms2_global.mod

ms2_stopwatch.mod: ms2_global.mod
