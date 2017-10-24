PROG =  read_data_file_dpaf
CSRCS = read_data_file_dpaf.c

OBJS =  $(CSRCS:.c=.o)	

NETINC = ../include
LIBDIR =
LIBS = -lm
CFLAGS =  -g -I$(NETINC)


debug := CFLAGS = -g

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(CFLAGS) $(OBJS) $(LIBDIR) $(LIBS) 

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

clean:
	-rm -f $(OBJS) tags core *.c.*

install:
	mv $(PROG) ../../bin/$(PROG)


