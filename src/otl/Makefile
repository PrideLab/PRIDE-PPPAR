# Makefile created by Yingda Deng (yingdadeng@whu.edu.cn), on jun 12, 2024


all: otl

DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ./
OBJ_PATH = .obj/

# Flags
FFLAGS = -g -O3
FTN = gfortran

LIBCOM = ../lib/libcom.a

SRC := $(wildcard *.f90)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))

otl : $(OBJ) $(LIBCOM)
	$(FTN) $(FFLAGS) $(OBJ) $(LIBCOM) -o $@

$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@

.PHONY : clean
clean :
	 -rm -rf otl $(OBJ_PATH)* $(OBJ_PATH)
