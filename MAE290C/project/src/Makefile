
FC = gfortran
FFTROOT = /usr/local
INCLUDE = $(FFTROOT)/include
LIB = $(FFTROOT)/lib
FLAGS = -lfftw3 -O3

OBJECTS = fft_mod.o \
          test_fft.o

MODULES = fft_mod.mod

.PHONY: clean

test_fft.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o test_fft.exe -I$(INCLUDE) -L$(LIB) $(FLAGS)

%.o: %.f90
	gfortran -c $<
%.mod: %.f90
	gfortran -c $< -I$(INCLUDE) -L$(LIB) $(FLAGS)

clean:
	rm -f $(OBJECTS) $(MODULES) test_fft.exe
