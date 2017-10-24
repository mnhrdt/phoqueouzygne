PROG =  ALOS_pre_process_SLC
CSRCS = ALOS_pre_process_SLC.c \
	parse_ALOS_commands.c \
	read_ALOS_data_SLC.c \
	swap_ALOS_data_info.c 


FSRCS =  
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../include

CLIBS = -L../lib -lALOS -lm
FLIBS = 
CFLAGS = -O3 
CFLAGS = -g -Wall 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)  $(FLIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
