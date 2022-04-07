
all: spp

DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ./
OBJ_PATH = .obj/

FFLAGS = -g -O3
FTN = gfortran

SRC := $(wildcard *.f90)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))

spp : $(OBJ)
	$(FTN) $(FFLAGS) $(OBJ) -o $@

$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@

.PHONY : clean
clean :
	-rm -rf spp $(OBJ_PATH)*


