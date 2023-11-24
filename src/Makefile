
SUBDIR = lib spp orbit tedit lsq redig arsig utils

all : $(SUBDIR)

$(SUBDIR) :
	make -j4 -C $@

install :
	@mkdir -p ../bin/
	mv arsig/arsig lsq/lsq tedit/tedit redig/redig orbit/sp3orb spp/spp \
		$(patsubst %.f90,%,$(wildcard utils/*.f90)) ../bin/

.PHONY : all $(SUBDIR) install
clean :
	-rm -rf $(addsuffix /.obj, $(SUBDIR))
	-rm -rf spp/*.o spp/*.mod
	-rm -rf lib/libcom.a
