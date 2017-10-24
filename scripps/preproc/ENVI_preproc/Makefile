LIB=./lib_src
ASA_CAT=./ASA_CAT
Dop_orbit=./Dop_orbit
ENVI_baseline=./ENVI_baseline
ENVI_decode=./ENVI_decode
ENVI_llt2rat=./ENVI_llt2rat
ENVI_look=./ENVI_look
SCRIPT=../ENVI_preproc/scripts

all:
	rm -r -f lib ../bin/ENV* ../bin/*env* ../bin/calc_dop_orb_envi
	mkdir lib 
	ln -s -f $(SCRIPT)/dump_orbit_envi.pl ../bin/.
	ln -s -f $(SCRIPT)/dump_time_envi.pl ../bin/.
	ln -s -f $(SCRIPT)/find_auxi.pl ../bin/.
	ln -s -f $(SCRIPT)/ENVI_pre_process ../bin/.
	cd $(LIB); make clean; make; make install; make clean
	cd $(ASA_CAT); make clean; make; make install; make clean
	cd $(Dop_orbit); make clean; make; make install; make clean
	cd $(ENVI_baseline); make clean; make; make install; make clean
	cd $(ENVI_decode); make clean; make; make install; make clean
	cd $(ENVI_llt2rat); make clean; make; make install; make clean
	cd $(ENVI_look); make clean; make; make install; make clean
