LIBDIR=../lib
CXXFLAGS= -O3 -std=c++11 -I. -fcx-fortran-rules -fcx-limited-range

target=$(LIBDIR)/libmodel_TopEffTh.a

all: $(target)

objects=HelAmps_TopEffTh.o Parameters_TopEffTh.o read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

