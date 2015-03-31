CXXFLAGS := -O3 -std=c++11 -Wall -l $(shell root-config --cflags) $(shell lhapdf-config --cflags)
LDFLAGS := -lm $(shell root-config --libs --glibs) -lFoam $(shell lhapdf-config --ldflags) -lcuba -lDelphes
CXX := g++

sourcedir := examples/

all: ttbar

ttbar: PROC_DIR := /home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum/
ttbar: MGproc := $(PROC_DIR)/SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.o $(PROC_DIR)/src/*.o

ttbar: $(sourcedir)/ME_ttbar_CUBA.o makefile
	$(CXX) $(LDFLAGS) $< $(MGproc) -o $(sourcedir)/ME_ttbar_CUBA.exe

$(sourcedir)/ME_ttbar_CUBA.o: $(sourcedir)/ME_ttbar_CUBA.cpp makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(PROC_DIR) -c $< -o $@

.PHONY: clean
clean:
	-rm $(sourcedir)/*.exe $(sourcedir)/*.o
