#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "fail.c"
#include "xmalloc.c"
#include "endianness.c"

#include "s1a.h"

struct bitstream {
	int byte;
	int bit;
	int total_bytes;
	uint8_t *data;
};

static void bitstream_init(struct bitstream *b, void *data, int total_bytes)
{
	b->byte = 0;
	b->bit = 0;
	b->data = data;
	b->total_bytes = total_bytes;
}

static bool bitstream_pop(struct bitstream *b)
{
	bool r = b->data[b->byte] & (1 << b->bit);

	b->byte += 1;
	b->bit = (b->bit + 1) % 8;
	if (b->byte >= b->total_bytes)
		fail("got too much bits from a bitstream!");

	return r;
}

void s1a_isp_sanity_check(struct s1a_isp *x)
{
	if (x->secondary_header.field.sync_marker != 0x352ef853)
		fail("incorrect sync marker found");

	if (x->secondary_header.field.baq_block_length != 31)
		fail("incorrect baq_block_field");
}

void s1a_decode_line(complex float *out, struct s1a_isp *x)
{
	s1a_isp_sanity_check(x);

	printf("S1A ISP:\n");
	printf("\tversion_number   = %d\n", x->version_number);
	printf("\tid               = %d\n", x->id);
	printf("\tsequence_control = %d\n", x->sequence_control);
	printf("\tdata_length      = %d\n", x->packet_data_length);

	printf("\tcoarse_T = %d\n", x->secondary_header.field.coarse_time);
	printf("\tfine_T   = %d\n", x->secondary_header.field.fine_time);
}
