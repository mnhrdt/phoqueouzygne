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

// Note: the following tables encode the Huffman trees on pages 71-73
// of the pdf named "S1-IF-ASD-PL-007" from the official Sentinel 1 docs.
//
// The trees are encoded by giving the state machine that traverses them.
// This form allows a very simple implementaiton of the decoder.
// States are numbered starting from 1, and in a breadht-first lexicographic
// ordering.  Leafs are indicated by numbers s<1 , and the symbol represented
// by a leaf is the value "-s".

static int huffman_table_brc_0[4][2] = {
	[1] = { 0,  2},
	[2] = {-1,  3},
	[3] = {-2, -3}
};
//       1 --- 2 --- 3 --- -3
//       |     |     |
//       |     |     |
//       0     -1    -2

static int huffman_table_brc_1[5][2] = {
	[1] =  { 0,  2},
	[2] =  {-1,  3},
	[3] =  {-2,  4},
	[4] =  {-3, -4}
};
//       1 --- 2 --- 3 --- 4 --- -4
//       |     |     |     |
//       |     |     |     |
//       0     -1    -2    -3

static int huffman_table_brc_2[7][2] = {
	[1] = { 0,  2},
	[2] = {-1,  3},
	[3] = {-2,  4},
	[4] = {-3,  5},
	[5] = {-4,  6},
	[6] = {-5, -6}
};
//       1 --- 2 --- 3 --- 4 --- 5 --- 6 --- -6
//       |     |     |     |     |     |
//       |     |     |     |     |     |
//       0     -1    -2    -3    -4    -5

static int huffman_table_brc_3[10][2] = {
	[1] = { 2,  3},
	[2] = { 0, -1},
	[3] = {-2,  4},
	[4] = {-3,  5},
	[5] = {-4,  6},
	[6] = {-5,  7},
	[7] = {-6,  8},
	[8] = {-7,  9},
	[9] = {-8, -9}
};
//       1 --- 3 --- 4 --- 5 --- 6 --- 7 --- 8 --- 9 --- -9
//       |     |     |     |     |     |     |     |
//       |     |     |     |     |     |     |     |
//       2     -2    -3    -4    -5    -6    -7    -8
//       |\
//       | \
//       0  -1

static int huffman_table_brc_4[16][2] = {
	[ 1] = {  2,   3},
	[ 2] = {  0,   4},
	[ 3] = {  5,   6},
	[ 4] = { -1,  -2},
	[ 5] = { -3,  -4},
	[ 6] = {  7,   8},
	[ 7] = { -5,  -6},
	[ 8] = { -7,   9},
	[ 9] = { -8,  10},
	[10] = { -9,  11},
	[11] = { 12,  13},
	[12] = {-10, -11},
	[13] = { 14,  15},
	[14] = {-12, -13},
	[15] = {-14, -15}
};
// (too difficult to draw in ascii)

static int (*huffman_table[5])[2] = {
	huffman_table_brc_0,
	huffman_table_brc_1,
	huffman_table_brc_2,
	huffman_table_brc_3,
	huffman_table_brc_4
};

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
