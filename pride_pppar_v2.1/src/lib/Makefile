# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn), on Nov 25, 2017
# Modified by Songfeng Yang (sfyang@whu.edu.cn) --- add multi-GNSS PPP

LIB := libcom.a
all : $(LIB)

DIR_GUARD = @mkdir -p $(@D)
OBJ_PATH = .obj/

# Flags
FTN = gfortran
FFLAGS = -g -O3 # -mcmodel=medium

SRC := $(wildcard *.f90)
SRC_FOR :=$(wildcard *.for)
SRC_F :=$(wildcard *.F)
OBJ := $(SRC:.f90=.o)
OBJ += $(SRC_FOR:.for=.o)
OBJ += $(SRC_F:.F=.o) 
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))
 
$(LIB) : $(OBJ)
	ar -rv $(LIB) $(OBJ)

$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@/bin/echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)

$(OBJ_PATH)%.o : %.for
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@/bin/echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)

$(OBJ_PATH)%.o : %.F
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@/bin/echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)
.PHONY : clean
clean :
	-rm $(LIB) $(OBJ_PATH)*

-include $(OBJ:.o=.d)

