########################################
# Makefile for ms2
########################################
MPI=1
ifneq ($(MPI),)
  F90 = mpif90
  CPPFLAGS = -DARCH=2 -DFORTRAN=90 -DMPI_VER=1 -D_PGF
  PROG =  ms2mpi
else
  F90 = pgf90 
  CPPFLAGS = -DARCH=2 -DFORTRAN=90 -DMPI_VER=0 -D_PGF
  PROG =  ms2
endif
F90FLAGS = -g -C
LDFLAGS =

SRCS =  ms2.F90 \
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

ms2_simulation.mod: ms2_ensemble.mod ms2_accumulator.mod ms2_global.mod

ms2_ensemble.mod: ms2_component.mod ms2_interaction.mod \
  ms2_accumulator.mod ms2_global.mod

ms2_component.mod: ms2_molecule.mod ms2_site.mod \
  ms2_accumulator.mod ms2_global.mod

ms2_molecule.mod: ms2_site.mod ms2_global.mod

ms2_site.mod: ms2_global.mod

ms2_interaction.mod: ms2_potential.mod ms2_component.mod ms2_global.mod

ms2_potential.mod: ms2_site.mod ms2_global.mod
