PROG =  offset_topo
CSRCS = offset_topo.c 

NETCDF = $(NETCDFHOME)/lib
NETINC = $(NETCDFHOME)/include
GMTLIB = $(GMTHOME)/lib
GMTINC = $(GMTHOME)/include
INC = ../../include 
FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)

CLIBS = -L../../lib -lgmtsar -lm
FLIBS = 
CFLAGS = -O2 -Bstatic $(M64)

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS) -L$(GMTLIB) -lgmt -L$(NETCDF) -lnetcdf

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC) -I$(GMTINC) -I$(NETINC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
