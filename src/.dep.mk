dump_headers.o: dump_headers.c xfopen.c fail.c xmalloc.c s1a.h
endianness.o: endianness.c
fail.o: fail.c
s1a_io.o: s1a_io.c fail.c xmalloc.c xfopen.c endianness.c s1a.h
xfopen.o: xfopen.c fail.c
xmalloc.o: xmalloc.c fail.c
