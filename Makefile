# C++ build
CXX = g++
CXXFLAGS = -std=c++11 -O2
CXXOBJS = gaussquad.o main.o
CXXTARGET = gaussquad

# Fortran build
FC = gfortran
FFLAGS = -O2 
FSRC = algama.f gaussquad_module.f90
F90SRC = ftest.f90
FTARGET = ftest

all: $(CXXTARGET) $(FTARGET)

$(CXXTARGET): $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(FTARGET): $(FSRC) $(F90SRC)
	$(FC) $(FFLAGS) $(FSRC) $(F90SRC) -o $(FTARGET)

clean:
	rm -f *.o *.mod $(CXXTARGET) $(FTARGET)
