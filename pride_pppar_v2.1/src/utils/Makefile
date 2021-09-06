# Modified by Songfeng Yang (sfyang@whu.edu.cn) --- add multi-GNSS PPP

all: xyz2enu enucov get_ctrl

DIR_GUARD = @mkdir -p $(@D)
OBJ_PATH = .obj/

# ASSIGMENTS
FTN = gfortran 
FFLAGS = -g -O3

LIBCOM = ../lib/libcom.a

xyz2enu : $(LIBCOM) xyz2enu.f90
	$(FTN) $(FFLAGS) xyz2enu.f90 $(LIBCOM) -o $@

enucov : $(LIBCOM) enucov.f90
	$(FTN) $(FFLAGS) enucov.f90 $(LIBCOM) -o $@

get_ctrl : $(LIBCOM) get_ctrl.f90
	$(FTN) $(FFLAGS) get_ctrl.f90 $(LIBCOM) -o $@


.PHONY : clean
clean :
	-rm xyz2enu enucov get_ctrl

