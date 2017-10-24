PROG =  read_sarleader_dpaf
CSRCS = read_sarleader_dpaf.c make_prm_dpaf.c

OBJS =  $(CSRCS:.c=.o)	

LIBS = -lm
SLIB = 
CFLAGS = -O3 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(LIBS) $(SLIB) 

.c.o:
	$(CC) -c $*.c -o $*.o

clean:
	-rm -f $(OBJS) tags core

install:
	mv $(PROG) ../../bin/$(PROG)

