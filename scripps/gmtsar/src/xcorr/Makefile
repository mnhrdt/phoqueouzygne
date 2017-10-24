PROG =  xcorr
CSRCS = xcorr.c \
	print_results.c \
	do_freq_xcorr.c \
	highres_corr.c \
	fft_interpolate_routines.c \
	parse_xcorr_input.c \
	get_locations.c \
	read_xcorr_data.c \
	cfft2d_1d_fftpack.c \
	do_time_int_xcorr.c

FSRCS =  
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../include

CLIBS = -L../../lib -lgmtsar -lm
FLIBS = 
CFLAGS = -O2 -Bstatic $(M64) -w

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS)  $(FLIBS) $(VEC)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
