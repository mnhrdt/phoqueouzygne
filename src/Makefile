CFLAGS  = -march=native -O3
LDLIBS  = -lm

OBJ     = s1a_io.o s1a_decode.o
BIN     = dump_headers decode_line

default : $(BIN)

$(BIN)  : $(OBJ)

clean   : ; $(RM) $(BIN) *.o

.dep.mk : ; $(CC) -MM *.c > $@
-include .dep.mk
