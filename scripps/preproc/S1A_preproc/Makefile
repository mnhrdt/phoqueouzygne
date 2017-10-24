#
#   make all the files
#
SRC=./src
LIBSAR=./libsar
TINYXPATH=./libtinyxpath

all:
	cd $(TINYXPATH); make clean; make
	cd $(LIBSAR); make clean; make
	cd $(SRC); make clean; make; make clean; make install
	cd $(TINYXPATH); make clean;
	cd $(LIBSAR); make clean;

