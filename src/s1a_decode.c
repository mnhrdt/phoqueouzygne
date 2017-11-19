#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "fail.c"
#include "xmalloc.c"
#include "endianness.c"

#include "s1a.h"

// data structure for bitwise memory traversal
struct bitstream {
	int byte;
	int bit;
	int total_bytes;
	uint8_t *data;
};

static void bitstream_init(struct bitstream *s, void *data, int total_bytes)
{
	s->byte = 0;
	s->bit = 0;
	s->data = data;
	s->total_bytes = total_bytes;
}

#include "smapa.h"
SMART_PARAMETER(REVERSE_BYTES,0)
SMART_PARAMETER(REVERSE_BITS,0)

static bool bitstream_pop(struct bitstream *s)
{

	int real_byte = s->byte;
	int real_bit = s->bit;

	if (REVERSE_BYTES()) {
		if (0 == s->byte%2)
			real_byte += 1;
		else
			real_byte -= 1;
	}

	if (REVERSE_BITS())
		real_bit = 7 - real_bit;

	assert(real_bit >= 0);
	assert(real_bit < 8);
	assert(real_byte >= 0);
	assert(real_byte < s->total_bytes);

	bool r = s->data[real_byte] & (1 << real_bit);
	s->byte += 1;
	s->bit = (s->bit + 1) % 8;
	if (s->byte >= s->total_bytes)
		fail("got too much bits from a bitstream!");

	printf("\tbit %d\n", r);
	return r;
}

static void bitstream_align_16(struct bitstream *s)
{
	int t[100], cx = 0;
	while (s->bit || s->byte % 2)
		t[cx++] = bitstream_pop(s);
	printf("%d lost bits:", cx);
	for (int i = 0; i < 0; i++)
		printf(" %d", t[i]);
	printf("\n");
}

// Note: the following tables encode the Huffman trees on pages 71-73
// of the pdf named "S1-IF-ASD-PL-007" from the official Sentinel 1 docs.
//
// The trees are encoded by giving the state machine that traverses them.
// This form allows a very simple implementaiton of the decoder.
// States are numbered starting from 1, and in a breadht-first lexicographic
// ordering.  Leafs are indicated by numbers s<1 , and the symbol represented
// by a leaf is the value "-s".

static int global_huffman_table_brc_0[4][2] = {
	[1] = { 0,  2},
	[2] = {-1,  3},
	[3] = {-2, -3}
};
//       1 --- 2 --- 3 --- -3
//       |     |     |
//       |     |     |
//       0     -1    -2

static int global_huffman_table_brc_1[5][2] = {
	[1] =  { 0,  2},
	[2] =  {-1,  3},
	[3] =  {-2,  4},
	[4] =  {-3, -4}
};
//       1 --- 2 --- 3 --- 4 --- -4
//       |     |     |     |
//       |     |     |     |
//       0     -1    -2    -3

static int global_huffman_table_brc_2[7][2] = {
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

static int global_huffman_table_brc_3[10][2] = {
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

static int global_huffman_table_brc_4[16][2] = {
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

static int (*global_huffman_table[5])[2] = {
	global_huffman_table_brc_0,
	global_huffman_table_brc_1,
	global_huffman_table_brc_2,
	global_huffman_table_brc_3,
	global_huffman_table_brc_4
};

void huffman_decode(struct bitstream *x)
{
	int brc = 3;
	int (*h)[2] = global_huffman_table[brc];
	int s = 1;
	while (1)
		if (s > 0)
			s = h[s][bitstream_pop(x)];
		else {
			printf("got symbol %d\n", -s);
	//		bytestream_push(y, -s);
			s = 1;
		}

}

void s1a_isp_verify_sanity(struct s1a_isp *x)
{
	if (x->secondary_header.field.sync_marker != 0x352ef853)
		fail("incorrect sync marker found");

	if (x->secondary_header.field.baq_block_length != 31)
		fail("incorrect baq_block_field");
}

void s1a_isp_verify_nominal_mode(struct s1a_isp *x)
{
	int TSTMOD = x->secondary_header.field.test_mode;
	int SIGTYP = x->secondary_header.field.signal_type;
	int BAQMOD = x->secondary_header.field.baq_mode;
	if (TSTMOD !=  0) fail("bad TSTMOD %d (should be 0)\n", TSTMOD);
	if (SIGTYP !=  0) fail("bad SIGTYP %d (should be 0)\n", SIGTYP);
	if (BAQMOD != 12) fail("bad BAQMOD %d (should be 12)\n", BAQMOD);
}

static void s1a_print_some_isp_fields(struct s1a_isp *x)
{
	printf("S1A ISP:\n");
	printf("\tversion_number   = %d\n", x->version_number);
	printf("\tid               = %d\n", x->id);
	printf("\tsequence_control = %d\n", x->sequence_control);
	printf("\tdata_length      = %d\n", x->packet_data_length);
	printf("\tdata_size        = %d\n", x->data_size);

	printf("\tcoarse_T = %d\n", x->secondary_header.field.coarse_time);
	printf("\tfine_T   = %d\n", x->secondary_header.field.fine_time);
	printf("\tDTID   = %d\n", x->secondary_header.field.data_take_id);
	printf("\tECC    = %d\n", x->secondary_header.field.ecc_number);
	printf("\tTSTMOD = %d\n", x->secondary_header.field.test_mode);
	printf("\tRXCHID = %d\n", x->secondary_header.field.rx_channel_id);
	printf("\tBAQMOD = %d\n", x->secondary_header.field.baq_mode);
	//printf("\tBAQLEN  = %d\n",x->secondary_header.field.baq_block_length);
	printf("\tRGDEC  = %d\n", x->secondary_header.field.range_decimation);
	printf("\tRXG    = %d\n", x->secondary_header.field.rx_gain);
	printf("\tTXPRR  = %d\n", x->secondary_header.field.tx_ramp_rate);
	printf("\tTXPSF  = %d\n", x->secondary_header.field.tx_pulse_start_frequency);
	printf("\tTXPL   = %d\n", x->secondary_header.field.tx_pulse_length);
	printf("\tRANK   = %d\n", x->secondary_header.field.rank);
	printf("\tPOL    = %d\n", x->secondary_header.field.polarisation);
	printf("\tCALMOD = %d\n", x->secondary_header.field.cal_mode);
	printf("\tTXPNO  = %d\n", x->secondary_header.field.tx_pulse_number);
	printf("\tSIGTYP = %d\n", x->secondary_header.field.signal_type);
	printf("\tSWAP   = %d\n", x->secondary_header.field.swap);
	printf("\tSWATH  = %d\n", x->secondary_header.field.swath_number);
	printf("\tNQ     = %d\n", x->secondary_header.field.number_of_quads);
}

void s1a_decode_line(complex float *out, struct s1a_isp *x)
{
	s1a_isp_verify_sanity(x);
	s1a_print_some_isp_fields(x);
	s1a_isp_verify_nominal_mode(x);

	// setup init bitstream
	struct bitstream s[1];
	bitstream_init(s, x->data, x->data_size);

	// space for output quads
	int NQ = x->secondary_header.field.number_of_quads;
	double out_ie[NQ];
	double out_io[NQ];
	double out_qe[NQ];
	double out_qo[NQ];

	// variables (cf. page 60 of document S1-IF-ASD-PL-0007)
	int NB = ceil(NQ / 128.0); // number of blocks
	int BRC[NB];               // bit rate of each block
	int code_ie[NQ];
	int code_io[NQ];
	int code_qe[NQ];
	int code_qo[NQ];
	//int b = 0;                 // block index (iterator)

	// decode IE data
	printf("NB = %d\n", NB);
	for (int b = 0; b < NB; b++)
	{
		BRC[b] = 4 * bitstream_pop(s)
		       + 2 * bitstream_pop(s)
		       + 1 * bitstream_pop(s);
		printf("BRC[%d] = %d\n", b, BRC[b]);
		if (BRC[b] < 0 || BRC[b] > 4) fail("bad BRC=%d\n", BRC[b]);
		int (*huf)[2] = global_huffman_table[BRC[b]];

		int num_hcodes = b < NB-1 ? 128 : NQ - 128 * (NB-1); // page 70
		printf("num_hcodes(%d) = %d\n", b, num_hcodes);
		for (int i = 0; i < num_hcodes; i++) // flowchart page 70
		{
			int sign = bitstream_pop(s) ? -1 : 1;
			int state = 1;
			while (state > 0)
				state = huf[state][bitstream_pop(s)];
			code_ie[i] = sign * (-state);
			printf("ie[%d] = %d\n", i, code_ie[i]);
		}
	}
	bitstream_align_16(s);
}
