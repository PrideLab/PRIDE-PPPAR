# Modified by Songfeng Yang (sfyang@whu.edu.cn) --- add multi-GNSS PPP

all: xyz2enu get_ctrl pbopos

DIR_GUARD = @mkdir -p $(@D)
OBJ_PATH = .obj/

# ASSIGMENTS
FTN = gfortran 
FFLAGS = -g -O3

LIBCOM = ../lib/libcom.a

xyz2enu : $(LIBCOM) xyz2enu.f90
	$(FTN) $(FFLAGS) xyz2enu.f90 $(LIBCOM) -o $@

pbopos : $(LIBCOM) pbopos.f90
	$(FTN) $(FFLAGS) pbopos.f90 $(LIBCOM) -o $@

get_ctrl : $(LIBCOM) get_ctrl.f90
	$(FTN) $(FFLAGS) get_ctrl.f90 $(LIBCOM) -o $@


.PHONY : clean
clean :
	-rm xyz2enu pbopos get_ctrl

