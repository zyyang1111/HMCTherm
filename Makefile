##################################################################################
#  HMCTemp v1.0 - 2017.06.26
#  A Cycle-accurate Simulator for Hybrid Memory Cube with Thermal analysis
#
#  Copyright (c) 2017, Zhiyuan Yang
#                      Ankur Srivastava
#                      University of Maryland 
#                      zyyang [at] umd [dot] edu
#  All rights reserved.
##################################################################################
PLAT = _PTHREAD

SuperLUroot	= $(HOME)/SuperLU_MT_3.1
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_mt$(PLAT).a
BLASLIB = $(SuperLUroot)/lib/libblas$(PLAT).a
LIBS		= $(SUPERLULIB) $(BLASLIB)
# HEADER = $(SuperLUroot)/SRC

CXXFLAGS=-O3 -std=c++11 -g -DDEBUG_LOG
EXE_NAME=HMCTherm
SRCDIR=sources

SRC = $(wildcard $(SRCDIR)/*.cpp)

OBJ = $(subst LibraryStubs.o,,$(addsuffix .o, $(basename $(SRC))))

REBUILDABLES=$(OBJ) $(EXE_NAME)

all: ${EXE_NAME}

#   $@ target name, $^ target deps, $< matched pattern
$(EXE_NAME): $(OBJ) $(SRCDIR)/thermal_solver.o
	$(CXX) $(LINK_FLAGS) -o $@ $^ $(LIBS) -lpthread -lm
	@echo "Built $@ successfully" 

#include the autogenerated dependency files for each .o file
-include $(OBJ:.o=.dep)

# build dependency list via gcc -M and save to a .dep file
%.dep : %.cpp
	@$(CXX) -M $(CXXFLAGS) $< > $@

# build all .cpp files to .o files
%.o : %.cpp
	g++ $(CXXFLAGS) -o $@ -c $<

$(SRCDIR)/thermal_solver.o: $(SRCDIR)/thermal_solver.c
	gcc -D__PTHREAD -g -D_LONGINT -DAdd_ -I/$(SuperLUroot)/SRC -o $@ -c $< 

$(SRCDIR)/sp_ienv.o: $(SRCDIR)/sp_ienv.c
	gcc -D__PTHREAD -g -D_LONGINT -DAdd_ -I/$(SuperLUroot)/SRC -o $@ -c $< 

clean: 
	-rm -f ${REBUILDABLES} $(SRCDIR)/*.dep ./result/*.log *.out *.csv ./power_trace/*.csv ./temperature_trace/*.csv $(SRCDIR)/*.o

remove:
	-rm -f ./result/*.log *.out *.csv ./power_trace/*.csv ./temperature_trace/*.csv
