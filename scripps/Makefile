#
#   make all the pre processors and gmtsar
#
ALOSPRE=./preproc/ALOS_preproc
ERSPRE=./preproc/ERS_preproc
ENVIPRE=./preproc/ENVI_preproc
S1APRE=./preproc/S1A_preproc
CSKPRE=./preproc/CSK_preproc
TSXPRE=./preproc/TSX_preproc
RS2PRE=./preproc/RS2_preproc
GMTSAR=./gmtsar
SNAPHU=./snaphu/src

all:
	cd $(ALOSPRE); make 
	cd $(ENVIPRE); make
	cd $(ERSPRE); make
	cd $(GMTSAR); make 
	cd $(SNAPHU); make clean ; make ; make clean
	cd $(S1APRE); make
#	cd $(TSXPRE); make
#	cd $(CSKPRE); make
#	cd $(RS2PRE); make
