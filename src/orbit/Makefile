# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn), on Dec 23, 2017
# Modified by Songfeng Yang (sfyang@whu.edu.cn) --- add multi-GNSS PPP

all: sp3orb

DIR_GUARD = @mkdir -p $(@D)
BIN_PATH = ./
OBJ_PATH = .obj/

# Flags
FFLAGS = -g -O3
FTN = gfortran

LIB = ../lib/libcom.a

sp3orb : $(LIB) sp3orb.f90 get_sp3orb_args.f90 lagrange_interp_sp3.f90 ../header/const.h ../header/orbit.h
	$(FTN) $(FFLAGS) sp3orb.f90 get_sp3orb_args.f90 lagrange_interp_sp3.f90 $(LIB) -o $@

.PHONY : clean
clean :
	-rm sp3orb

