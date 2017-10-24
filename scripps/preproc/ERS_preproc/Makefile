ers_line_fixer=./ers_line_fixer
read_data_file_dpaf=./read_data_file_dpaf
read_data_file_ccrs=./read_data_file_ccrs
read_sarleader_dpaf=./read_sarleader_dpaf
SCRIPT=../ERS_preproc/scripts

all:
	rm -r -f ../bin/ERS* ../bin/*ers* ../bin/read*
	ln -s -f $(SCRIPT)/dump_orbit_ers.pl ../bin/.
	ln -s -f $(SCRIPT)/ERS_pre_process ../bin/.
	ln -s -f ../bin/ENVI_baseline ../bin/ERS_baseline
	ln -s -f ../bin/ENVI_llt2rat  ../bin/ERS_llt2rat
	cd $(ers_line_fixer); make clean; make; make install; make clean
	cd $(read_data_file_dpaf); make clean; make; make install; make clean
	cd $(read_data_file_ccrs); make clean; make; make install; make clean
	cd $(read_sarleader_dpaf); make clean; make; make install; make clean
