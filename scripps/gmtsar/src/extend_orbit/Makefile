PROG =  extend_orbit
CSRCS = extend_orbit.c 

FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../include

CLIBS = -L../../lib -lgmtsar -lm
FLIBS = 
CFLAGS = -O2 -Bstatic $(M64)

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
