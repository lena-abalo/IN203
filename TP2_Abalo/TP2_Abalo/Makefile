include Make.inc
#include Make_msys2.inc
#include Make_osx.inc

CXXFLAGS = -std=c++11 -fPIC
ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
else
CXXFLAGS += -O3 -march=native -Wall
endif

ALL=matvec_col.exe Mandelbrot.exe Mandelbrot_me.exe test.exe matvec_lig.exe

default:	help

all: $(ALL)

clean:
	@rm -fr *.o *.exe *~ *.tga

%.exe: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB)


help:
	@echo "Available targets : "
	@echo "    all            : compile all executables"
	@echo "    matvec_col.exe     : compile matrix vector col product executable"
	@echo "    matvec_lig.exe     : compile matrix vector lig product executable"
	@echo "    test.exe     : compile test executable"
	@echo "    Mandelbrot.exe : compile Mandelbrot set computation executable"
	@echo "    Mandelbrot_me.exe : compile Mandelbrot maitre esclave set computation executable"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"
