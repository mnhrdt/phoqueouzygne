# Makefile for lib_src 
LIB = libsar.a
INC = -I../include -I$(NETCDFHOME)/include
#INC = -I../include -I$(GMTHOME)/include -I$(NETCDFHOME)/include
CXX = g++ 
CXXFLAGS = -O2  -Bstatic $(INC) $(M64) -w

SRCS =	doppler.cpp \
	gmtsar.cpp \
	instrument.cpp \
	orbit.cpp \
	planet.cpp \
	platform.cpp \
	scene.cpp \
	statevector.cpp \
	utilities.cpp

OBJS= $(SRCS:.cpp=.o)

$(LIB) : $(OBJS)
	ar rv $(LIB) $?
	ranlib  $(LIB)

clean:
	\rm -f *.o 
