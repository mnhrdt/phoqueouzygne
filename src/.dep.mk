bits.o: bits.c
btt.o: btt.c iio.h xmalloc.c fail.c
decode_line.o: decode_line.c s1a.h
dump_headers.o: dump_headers.c xfopen.c fail.c xmalloc.c s1a.h
extract_img.o: extract_img.c s1a.h xmalloc.c fail.c xfopen.c smapa.h \
 iio.h pickopt.c
fail.o: fail.c
iio.o: iio.c
pickopt.o: pickopt.c
s1a_decode.o: s1a_decode.c fail.c bits.c s1a.h
s1a_focus.o: s1a_focus.c s1a.h smapa.h
s1a_io.o: s1a_io.c fail.c xmalloc.c xfopen.c bits.c s1a.h
xfopen.o: xfopen.c fail.c
xmalloc.o: xmalloc.c fail.c
