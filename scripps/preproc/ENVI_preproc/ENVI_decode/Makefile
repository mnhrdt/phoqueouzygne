PROG =  asa_im_decode
CSRCS = asa_im_decode.c 

FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = 

CLIBS = -lm
FLIBS = 
CFLAGS = -O3 -w

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
