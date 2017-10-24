PROG =  phase2topo
CSRCS = phase2topo.c 

NETCDF = $(NETCDFHOME)/lib
GMTLIB = $(GMTHOME)/lib
FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../include

CLIBS = -L../../lib -lgmtsar -lm
FLIBS = 
CFLAGS = -O2 -Bstatic $(M64) -w

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS) -L$(GMTLIB) -lgmt -L$(NETCDF) -lnetcdf

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
