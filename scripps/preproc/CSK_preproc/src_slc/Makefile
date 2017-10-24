PROG =  make_slc_csk
CSRCS = make_slc_csk.cpp csk.cpp

FSRCS =
FC = 
CC = g++

OBJS =  $(CSRCS:.cpp=.o)
INC = -I../include -I../libsar/include -I../include -I$(NETCDFHOME)/include

CLIBS = -L../libtinyxpath -ltinyxpath -L../libsar -lsar -L$(NETCDFHOME)/lib -lboost_system$(mt) -lboost_program_options$(mt) -lboost_filesystem$(mt) -lboost_date_time$(mt) -lhdf5_cpp -lhdf5 -lm
CXXFLAGS = -O2

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) 

make_slc_csk.o:
	$(CC) $(CXXFLAGS) -c $(CSRCS) $(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
