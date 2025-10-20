# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn), on Mar 9, 2018
# Modified by Songfeng Yang (sfyang@whu.edu.cn) --- add multi-GNSS PPP

all: arsig

DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ./
OBJ_PATH = .obj/

# ASSIGMENTS
FTN = gfortran
FFLAGS = -g -O0
CFLAGS = -std=c99

LIBCOM = ../lib/libcom.a

SRC := $(wildcard *.f90)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))
COBJ := $(OBJ_PATH)xgb_model.o

arsig : $(OBJ) $(COBJ) $(LIBCOM)
	$(FTN) $(FFLAGS) $(OBJ) $(COBJ) $(LIBCOM) -lm -o $@

$(LIB) : $(OBJ)
	ar -rv $(LIB) $(OBJ)
 
$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@/bin/echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)
	
$(OBJ_PATH)%.o : %.c
	$(DIR_GUARD)
	gcc $(CFLAGS) -c $< -o $@


 
.PHONY : clean
clean :
	-rm $(OBJ_PATH)*
 
 -include $(OBJ:.o=.d)

