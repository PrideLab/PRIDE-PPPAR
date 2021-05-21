#
OBJ	= spp.o rtklib.o rtkcmn.o rinex.o rtkpos.o postpos.o solution.o pntpos.o ephemeris.o

CFLAGS	= -O3 -g

all	: spp

spp	: $(OBJ)
	gfortran $(CFLAGS) $^ -o $@
	chmod +x spp

spp.o : spp.f90 postpos.o
	gfortran -c $(CFLAGS) $< -o $@
postpos.o : postpos.f90 rtkpos.o solution.o rinex.o
	gfortran -c $(CFLAGS) $< -o $@
rinex.o : rinex.f90 rtkcmn.o
	gfortran -c $(CFLAGS) $< -o $@
solution.o : solution.f90 rtkcmn.o
	gfortran -c $(CFLAGS) $< -o $@
rtkpos.o : rtkpos.f90 pntpos.o
	gfortran -c $(CFLAGS) $< -o $@
pntpos.o : pntpos.f90 ephemeris.o
	gfortran -c $(CFLAGS) $< -o $@
ephemeris.o : ephemeris.f90 rtkcmn.o
	gfortran -c $(CFLAGS) $< -o $@
rtkcmn.o : rtkcmn.f90 rtklib.o
	gfortran -c $(CFLAGS) $< -o $@
rtklib.o : rtklib.f90
	gfortran -c $(CFLAGS) $< -o $@

.PHONY:clean
clean:
	rm -rf *.o *.mod
