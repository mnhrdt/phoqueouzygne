# Makefile for lib src 
CC = gcc 
LIB = libENVI.a
INC = ../../ALOS_preproc/include
CFLAGS = -Wall -I$(INC)

SRCS =	read_ENVI_orb.c \
	ENVI_ldr_orbit.c

OBJS= $(SRCS:.c=.o)

$(LIB) : $(OBJS)
	ar r $(LIB) $?
	ranlib  $(LIB)

install:
	mv *.a ../lib
	ranlib  ../lib/$(LIB)

clean:
	\rm -f *.a *.o 
