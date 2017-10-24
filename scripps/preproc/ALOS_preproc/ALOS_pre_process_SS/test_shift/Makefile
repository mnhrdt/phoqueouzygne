PROG =  test_shift
CSRCS = test_shift.c \
	shift.c
OS = `arch`

FSRCS =  
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../include

CLIBS = -L../../lib/$(OS) -lALOS -lm
FLIBS = 
CFLAGS = -O3 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)  $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../bin/$(OS)
clean:
	-rm -f $(OBJS) tags core
