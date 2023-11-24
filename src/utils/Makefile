# Makefile created by Jihang Lin (jihanglin@whu.edu.cn), on Feb 23, 2023 

# Library (should be compiled first)
LIBCOM = ../lib/libcom.a

# Flags
FTN = gfortran 
FFLAGS = -g -O3

# Get all Fortran files in the current directory
F90_FILES = $(wildcard *.f90)
EXE_FILES = $(patsubst %.f90,%,$(F90_FILES))

all : $(EXE_FILES)

%: %.f90 $(LIBCOM) 
	$(FTN) $(FFLAGS) $< $(LIBCOM) -o $@

.PHONY : all clean
clean :
	rm -f $(EXE_FILES) 
