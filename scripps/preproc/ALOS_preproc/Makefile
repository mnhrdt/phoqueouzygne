LIB=./lib_src
ALOS_pre_process=./ALOS_pre_process
ALOS_baseline=./ALOS_baseline
ALOS_merge=./ALOS_merge
ALOS_fbd2fbs=./ALOS_fbd2fbs
ALOS_fbd2fbs_SLC=./ALOS_fbd2fbs_SLC
ALOS_fbd2ss=./ALOS_fbd2ss
ALOS_llt2rat=./ALOS_llt2rat
ALOS_pre_process_SS=./ALOS_pre_process_SS
ALOS_pre_process_SLC=./ALOS_pre_process_SLC
ALOS_look=./ALOS_look

all:
	rm -r -f lib ../bin/ALOS*
	mkdir lib 
	cd $(LIB); make clean; make; make install; make clean
	cd $(ALOS_pre_process); make clean; make; make install; make clean
	cd $(ALOS_baseline); make clean; make; make install; make clean
	cd $(ALOS_merge); make clean; make; make install; make clean
	cd $(ALOS_fbd2fbs); make clean; make; make install; make clean
	cd $(ALOS_fbd2fbs_SLC); make clean; make; make install; make clean
	cd $(ALOS_fbd2ss); make clean; make; make install; make clean
	cd $(ALOS_llt2rat); make clean; make; make install; make clean
	cd $(ALOS_pre_process_SS); make clean; make; make install; make clean
	cd $(ALOS_pre_process_SLC); make clean; make; make install; make clean
	cd $(ALOS_look); make clean; make; make install; make clean
