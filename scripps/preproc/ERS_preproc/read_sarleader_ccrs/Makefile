PROG =  read_sarleader_ccrs
CSRCS = read_sarleader_ccrs.c make_prm_ccrs.c

OBJS =  $(CSRCS:.c=.o)	

LIBS = -lm
SLIB = 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(LIBS) $(SLIB) 

.c.o:
	$(CC) -c $*.c -o $*.o

clean:
	-rm -f $(OBJS) tags core

install:
	mv $(PROG) ../../bin/$(PROG)

