CFLAGS  = -march=native -O3
LDLIBS  = -lfftw3f -lm

OBJ     = s1a_io.o s1a_decode.o s1a_focus.o
BIN     = dump_headers decode_line extract_img

default : $(BIN)

$(BIN)  : $(OBJ)

test    : $(BIN)
	./decode_line ../data/vv.dat 2000

clean   : ; $(RM) $(BIN) *.o

.dep.mk : ; $(CC) -MM *.c > $@
-include .dep.mk
