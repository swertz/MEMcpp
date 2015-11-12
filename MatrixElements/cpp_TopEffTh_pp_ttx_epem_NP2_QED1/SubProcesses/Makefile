LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_TopEffTh
CXXFLAGS= -std=c++11 -O3 -fcx-fortran-rules -fcx-limited-range
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

all: cpp_TopEffTh_pp_ttx_epem_NP2_QED1.o

cpp_TopEffTh_pp_ttx_epem_NP2_QED1.o: cpp_TopEffTh_pp_ttx_epem_NP2_QED1.cc $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LIBFLAGS) -I$(INCDIR)

$(LIBDIR)/lib$(MODELLIB).a:
	pushd ../../src && make && popd

.PHONY: clean

clean:
	rm -f *.o
	pushd ../../src && make clean && popd
