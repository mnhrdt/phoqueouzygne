CFLAGS  = -march=native -O3 -Wall -Wno-unused
LDLIBS  = -lfftw3f -lm -ljpeg -lpng -ltiff

OBJ     = s1a_io.o s1a_decode.o s1a_focus.o iio.o
BIN     = dump_headers decode_line extract_img unsubcommutate btt

default : $(BIN)

$(BIN)  : $(OBJ)

test    : $(BIN)
	./decode_line ../data/vv.dat 2000

clean   : ; $(RM) $(BIN) *.o

.dep.mk : ; $(CC) -MM *.c > $@
-include .dep.mk
