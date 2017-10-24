LIB=./lib_src
esarp=./src/esarp
conv=./src/conv
xcorr=./src/xcorr
phasediff=./src/phasediff
offset_topo=./src/offset_topo
phasefilt=./src/phasefilt
phase2topo=./src/phase2topo
bperp=./src/bperp
resamp=./src/resamp
calc_dop_orb=./src/calc_dop_orb
baseline=./src/baseline
llt2rat=./src/llt2rat
look=./src/look
extend_orbit=./src/extend_orbit
sbas=./src/sbas


all:
	rm -r -f lib bin
	mkdir lib bin
	cd $(LIB); make clean; make; make install
	cd $(esarp); make clean; make; make install; make clean
	cd $(conv); make clean; make; make install; make clean
	cd $(xcorr); make clean; make; make install; make clean
	cd $(phasediff); make clean; make; make install; make clean
	cd $(offset_topo); make clean; make; make install; make clean
	cd $(phasefilt); make clean; make; make install; make clean
	cd $(phase2topo); make clean; make; make install; make clean
	cd $(bperp); make clean; make; make install; make clean
	cd $(resamp); make clean; make; make install; make clean
	cd $(calc_dop_orb); make clean; make; make install; make clean
	cd $(baseline); make clean; make; make install; make clean
	cd $(llt2rat); make clean; make; make install; make clean
	cd $(look); make clean; make; make install; make clean
	cd $(extend_orbit); make clean; make; make install; make clean
#	cd $(sbas); make clean; make; make install; make clean
	cd $(LIB); make clean
