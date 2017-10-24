PROG =  calc_dop_orb_envi
CSRCS = calc_dop_orb_envi.c 

FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../ALOS_preproc/include

CLIBS = -L../../ALOS_preproc/lib -L../lib -lENVI -lALOS -lm
FLIBS = 
CFLAGS = -O3 -Wall 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
