CPP= g++
CC = gcc
CPPFLAGS=
CCFLAGS=
F90FLAGS=
F90 = gfortran

FRAMEWORK_INC=/home/hpcdjenz/workspace/Steereo/src
FRAMEWORK_LIB=/home/hpcdjenz/workspace/Steereo/src
FRAMEWORK_NAME=Steereo

LIBS= -lpthread -l$(FRAMEWORK_NAME) -L$(FRAMEWORK_LIB) -lc -lstdc++
INCLUDE= -I$(FRAMEWORK_INC)
MODULES= megaMolSnapshotCommand.so

all: $(MODULES)

%.so: %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -fPIC -shared -Wl,-soname,$@ -o $@ $^

clean:
	rm $(MODULES)