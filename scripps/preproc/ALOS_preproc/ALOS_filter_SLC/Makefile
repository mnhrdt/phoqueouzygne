PROG =  ALOS_filter_SLC
CSRCS = ALOS_filter_SLC.c

FSRCS =  
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../include

CLIBS = -L../lib -lALOS -lm
CFLAGS = -O3 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(VEC)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
