PROG =  resamp
CSRCS = resamp.c \
	print_fix_prm.c \
	ram2ras.c \
	interp.c \
	bicubic_one.c\
	cubic_kernel.c\
	sinc_one.c \
	sinc_kernel.c

NETCDF = $(NETCDFHOME)/lib
GMTLIB = $(GMTHOME)/lib
FSRCS = 
FC = gfortran
CC = gcc

OBJS =  $(CSRCS:.c=.o) $(FSRCS:.f=.o)	

INC = ../../include

CLIBS = -L../../lib -lgmtsar -lm
FLIBS =
CFLAGS = -Bstatic $(M64) -w
FFLAGS = -Bstatic $(M64)

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS) -L$(GMTLIB) -lgmt -L$(NETCDF) -lnetcdf

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

clean:
	-rm -f $(OBJS) tags core

install:
	mv $(PROG) ../../bin
