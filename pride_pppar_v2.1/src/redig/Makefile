# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn) on Mar 8 2018
# Modified by Songfeng Yang (sfyang@whu.edu.cn) --- add multi-GNSS PPP

all: redig

DIR_GUARD = @mkdir -p $(@D)
OBJ_PATH = .obj/

# ASSIGMENTS
FTN = gfortran 
FFLAGS = -g -O3

LIBCOM = ../lib/libcom.a  

SRC := $(wildcard *.f90)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))

redig : $(OBJ) $(LIBCOM)
	$(FTN) $(FFLAGS) $(OBJ) $(LIBCOM) -o $@

$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@/bin/echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)

.PHONY : clean
clean :
	-rm $(OBJ_PATH)*

-include $(OBJ:.o=.d)

