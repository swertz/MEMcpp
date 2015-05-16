CXXFLAGS := -std=c++11 -O2 -g -Wall $(shell root-config --cflags) $(shell lhapdf-config --cflags)
LDFLAGS := -lm $(shell root-config --libs --glibs) $(shell lhapdf-config --ldflags) -lcuba -lDelphes
CXX := g++

sourcedir := examples/

#### TTbar specific

ttbar_PROC_DIR := /home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum/
ttbar_MGproc := $(ttbar_PROC_DIR)/SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.o
#$(ttbar_PROC_DIR)/src/*.o


all: ttbar



#### TTbar specific

ttbar: $(sourcedir)/ME_ttbar_CUBA.o $(sourcedir)/utils.o $(sourcedir)/jacobianD.o $(sourcedir)/MEEvent.o $(sourcedir)/Integrand_TTbar.o $(sourcedir)/binnedTF.o $(sourcedir)/transferFunction.o $(sourcedir)/MEWeight.o $(ttbar_MGproc)
	$(CXX) -o $(sourcedir)/ME_ttbar_CUBA.exe $^ $(LDFLAGS) -L$(ttbar_PROC_DIR)/lib/ -lmodel_sm

$(sourcedir)/ME_ttbar_CUBA.o: $(sourcedir)/ME_ttbar_CUBA.cpp
	$(CXX) $(CXXFLAGS)  -I$(ttbar_PROC_DIR) -c $< -o $@

$(sourcedir)/MEEvent.o: $(sourcedir)/MEEvent.cpp $(sourcedir)/MEEvent.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(sourcedir)/MEWeight.o: $(sourcedir)/MEWeight.cpp $(sourcedir)/MEWeight.h
	$(CXX) $(CXXFLAGS) -I$(ttbar_PROC_DIR) -c $< -o $@

$(sourcedir)/Integrand_TTbar.o: $(sourcedir)/Integrand_TTbar.cpp $(sourcedir)/MEWeight.h
	$(CXX) $(CXXFLAGS) -I$(ttbar_PROC_DIR) -c $< -o $@

##### Common part

$(sourcedir)/jacobianD.o: $(sourcedir)/jacobianD.cpp $(sourcedir)/jacobianD.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(sourcedir)/utils.o: $(sourcedir)/utils.cpp $(sourcedir)/utils.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(sourcedir)/binnedTF.o: $(sourcedir)/binnedTF.cpp $(sourcedir)/binnedTF.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(sourcedir)/transferFunction.o: $(sourcedir)/transferFunction.cpp $(sourcedir)/transferFunction.h
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	-rm $(sourcedir)/*.exe $(sourcedir)/*.o
