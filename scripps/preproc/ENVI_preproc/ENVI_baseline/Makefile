PROG =  ENVI_baseline
CSRCS = ENVI_baseline.c \
	ENVI_llt2rat_sub.c

FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../ALOS_preproc/include

CLIBS = -L../../ALOS_preproc/lib -lALOS -L../lib -lENVI -lm
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
