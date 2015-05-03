
FC = gfortran
FFTROOT = /usr/local/Cellar/fftw/3.3.4_1
INCLUDE = $(FFTROOT)/include
LIB = $(FFTROOT)/lib
FLAG = -lfftw3

.PHONY = run 

pb3.exe : pb3.o
		$(FC) pb3.o -o $@ -L$(LIB) $(FLAG)

pb3.o : pb3.f90
		$(FC) -O -c pb3.f90 -I$(INCLUDE) -L$(LIB) $(FLAG)

run :
		./pb3



