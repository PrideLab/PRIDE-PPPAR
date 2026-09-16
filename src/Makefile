
SUBDIR = lib spp orbit tedit lsq redig arsig utils otl mhm

all : $(SUBDIR)

$(SUBDIR) :
	make -j4 -C $@

# Every other directory links ../lib/libcom.a, so a parallel top-level make
# (make -jN) must finish lib first.
$(filter-out lib,$(SUBDIR)) : lib

install :
	@mkdir -p ../bin/
	mv arsig/arsig lsq/lsq tedit/tedit redig/redig orbit/sp3orb spp/spp otl/otl mhm/mhm \
		$(patsubst %.f90,%,$(wildcard utils/*.f90)) ../bin/

.PHONY : all $(SUBDIR) install
clean :
	-rm -rf $(addsuffix /.obj, $(SUBDIR))
	-rm -rf spp/*.o spp/*.mod
	-rm -rf lib/libcom.a
	-rm -rf mhm/*.o
