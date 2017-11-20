decode_line.o: decode_line.c s1a.h
dump_headers.o: dump_headers.c xfopen.c fail.c xmalloc.c s1a.h
endianness.o: endianness.c
extract_img.o: extract_img.c s1a.h xmalloc.c fail.c xfopen.c
fail.o: fail.c
s1a_decode.o: s1a_decode.c fail.c xmalloc.c endianness.c s1a.h
s1a_io.o: s1a_io.c fail.c xmalloc.c xfopen.c endianness.c s1a.h
xfopen.o: xfopen.c fail.c
xmalloc.o: xmalloc.c fail.c
