PROG =  ALOS_fbd2ss
CSRCS = ALOS_fbd2ss.c

FSRCS =  
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../include

#CLIBS = -L../lib -lALOS -framework vecLib -lm
CLIBS = -L../lib -lALOS -lm
FLIBS = 
CFLAGS = -O3 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)  $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
