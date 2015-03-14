BOOST_DIR := /usr/include/boost141

DELPHES_DIR := /home/fynu/swertz/storage/Delphes/Delphes-3.1.2/

CXXFLAGS := -O3 -std=c++0x -Wall -l $(shell root-config --cflags) $(shell lhapdf-config --cflags) -I$(BOOST_DIR) -I$(DELPHES_DIR)
LDFLAGS := $(shell root-config --libs --glibs) -lFoam $(shell lhapdf-config --ldflags) -L$(DELPHES_DIR) -lDelphes
CXX := g++


sourcedir := examples/

all: ww ttbar

ttbar: PROC_DIR := /home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum/
ttbar: MGproc := $(PROC_DIR)/SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.o $(PROC_DIR)/src/*.o

ttbar: $(sourcedir)/ME_ttbar_test.o makefile
	$(CXX) $(LDFLAGS) $< $(MGproc) -o $(sourcedir)/ME_ttbar_test.exe

$(sourcedir)/ME_ttbar_test.o: $(sourcedir)/ME_ttbar_test.cpp makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(PROC_DIR) -c $< -o $@


ww: PROC_DIR := /home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/
ww: MGproc := $(PROC_DIR)/SubProcesses/P0_Sigma_sm_uux_epvemumvmx/CPPProcess.o $(PROC_DIR)/src/*.o

ww: $(sourcedir)/ME_ww_test.o makefile
	$(CXX) $(LDFLAGS) $< $(MGproc) -o $(sourcedir)/ME_ww_test.exe

$(sourcedir)/ME_ww_test.o: $(sourcedir)/ME_ww_test.cpp makefile
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(PROC_DIR) -c $< -o $@

#w: $(sourcedir)/ME_w_test.o makefile
#	export MGproc = $(PROC_DIR)/SubProcesses/P0_Sigma_sm_uux_epvemumvmx/CPPProcess.o $(PROC_DIR)/src/*.o
#	$(CXX) $(LDFLAGS) $< $(MGproc) -o $(sourcedir)/ME_w_test.exe
#
#$(sourcedir)/MW_w_test.o: $(sourcedir)/MW_w_test.cpp makefile
#	export PROC_DIR := /home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(PROC_DIR) -c $< -o $@


.PHONY: clean
clean:
	-rm $(sourcedir)/*.exe $(sourcedir)/*.o
