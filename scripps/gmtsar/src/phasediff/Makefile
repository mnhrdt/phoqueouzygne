PROG =  phasediff
CSRCS = phasediff.c \
	utils_complex.c \
	file_stuff.c \
	read_optional_args.c

NETCDF = $(NETCDFHOME)/lib
GMTLIB = $(GMTHOME)/lib
FSRCS = 
FC = gfortran
CC = gcc

OBJS =  $(CSRCS:.c=.o) $(FSRCS:.f=.o)	

INC = ../../include

CLIBS = -L../../lib -lgmtsar -lm
FLIBS =
CFLAGS = -O2 -Bstatic $(M64) -w
FFLAGS = -O2 -Bstatic $(M64)

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS) -L$(GMTLIB) -lgmt -L$(NETCDF) -lnetcdf

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

clean:
	-rm -f $(OBJS) tags core

install:
	mv $(PROG) ../../bin
