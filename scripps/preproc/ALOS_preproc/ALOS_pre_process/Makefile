PROG =  ALOS_pre_process
CSRCS = ALOS_pre_process.c \
	parse_ALOS_commands.c \
	read_ALOS_data.c \
	read_ALOSE_data.c \
	swap_ALOS_data_info.c \
	roi_utils.c

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
