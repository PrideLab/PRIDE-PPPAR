
SUBDIR = lib \
	 orbit \
	 lsq \
	 arsig \
	 tedit \
	 redig \
	 utils \
	 spp

all : $(SUBDIR)

$(SUBDIR) :
	make -j4 -C $@

install :
	@mkdir -p ../bin/
	mv arsig/arsig lsq/lsq tedit/tedit redig/redig orbit/sp3orb \
		utils/xyz2enu utils/enucov utils/get_ctrl spp/spp ../bin/

.PHONY : all $(SUBDIR) install

clean :
	-rm -rf $(addsuffix /.obj, $(SUBDIR))
	-rm -rf spp/*.o spp/*.mod
	-rm -rf lib/libcom.a

