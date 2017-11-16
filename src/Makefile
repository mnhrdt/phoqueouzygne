# configuration
CFLAGS  = -march=native -O3
LDLIBS  = -lm

# variables
OBJ     = s1a_io.o
BIN     = dump_headers

# default targets
default : $(BIN)

# dependences
$(BIN)  : $(OBJ)

# bureaucracy
clean   : ; $(RM) $(BIN) *.o

# automatic dependences (necessary for directly-included files .c files)
.dep.mk : ; $(CC) -MM *.c > $@
-include .dep.mk
