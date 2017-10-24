PROG =  make_slc_tsx
CSRCS = make_slc_tsx.cpp tsx.cpp Burst.cpp Cosar.cpp Header.cpp

FSRCS =
FC = 
CC = g++

OBJS =  $(CSRCS:.cpp=.o)
INC = -I../include -I../libsar/include -I../include -I$(NETCDFHOME)/include

CLIBS = -L../libtinyxpath -ltinyxpath -L../libsar -lsar -L$(NETCDFHOME)/lib -lboost_system$(mt) -lboost_program_options$(mt) -lboost_filesystem$(mt) -lboost_date_time$(mt) -ltiff -lm
CXXFLAGS = -O2 -w

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) 

make_slc_tsx.o:
	$(CC) $(CXXFLAGS) -c $(CSRCS) $(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
