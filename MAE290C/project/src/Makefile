
FC = gfortran
FFTROOT = /usr/local
INCLUDE = $(FFTROOT)/include
LIB = $(FFTROOT)/lib
FLAGS = -lfftw3_omp -lfftw3 -lm -fopenmp -O3

OBJECTS = fft_mod.o \
          vort_eqn.o

MODULES = fft_mod.mod

.PHONY: clean

#
# *** model ***
#

model.exe:	$(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o model.exe -I$(INCLUDE) -L$(LIB) $(FLAGS)

#
#  *** tests ***
#

test_wavenumber.exe:	$(MODULES) test_wavenumber.o
	gfortran $(OBJECTS) test_wavenumber.o -o test_wavenumber.exe -I$(INCLUDE) -L$(LIB) $(FLAGS)

test_fft.exe:	$(MODULES) test_fft.o
	gfortran $(OBJECTS) test_fft.o -o test_fft.exe -I$(INCLUDE) -L$(LIB) $(FLAGS)

%.o: %.f90
	gfortran -c $<
%.mod: %.f90
	gfortran -c $< -I$(INCLUDE) -L$(LIB) $(FLAGS)

clean:
	rm -f $(OBJECTS) $(MODULES) test_wavenumber.exe test_fft.exe model.exe
