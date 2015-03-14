BOOST_DIR := /usr/include/boost141

DELPHES_DIR=/home/fynu/swertz/storage/Delphes/Delphes-3.1.2/

PROC_DIR := /home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/

CXXFLAGS := -O3 -std=c++0x -Wall -l $(shell root-config --cflags) $(shell lhapdf-config --cflags) -I$(BOOST_DIR) -I$(DELPHES_DIR) -I$(PROC_DIR)
LDFLAGS := $(shell root-config --libs --glibs) -lFoam $(shell lhapdf-config --ldflags) -L$(DELPHES_DIR) -lDelphes
CXX := g++

MGproc := $(PROC_DIR)/SubProcesses/P0_Sigma_sm_uux_epvemumvmx/CPPProcess.o $(PROC_DIR)/src/*.o

sourcedir := examples/

all: ww w

ww: $(sourcedir)/ME_ww_test.o makefile
	$(CXX) $(LDFLAGS) $< $(MGproc) -o $(sourcedir)/ME_ww_test.exe

$(sourcedir)/MW_ww_test.o: $(sourcedir)/MW_ww_test.cpp makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@


w: $(sourcedir)/ME_w_test.o makefile
	$(CXX) $(LDFLAGS) $< $(MGproc) -o $(sourcedir)/ME_w_test.exe

$(sourcedir)/MW_w_test.o: $(sourcedir)/MW_w_test.cpp makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@


.PHONY: clean
clean:
	-rm $(sourcedir)/*.exe $(sourcedir)/*.o
