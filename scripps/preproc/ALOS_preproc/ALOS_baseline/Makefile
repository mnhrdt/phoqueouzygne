PROG =  ALOS_baseline
CSRCS = ALOS_baseline.c \
	ALOS_llt2rat_sub.c

FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../include

CLIBS = -L../lib -lALOS -lm
FLIBS = 
CFLAGS = -O3 -Wall -g

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
