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

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
		  (byte & 0x80 ? '1' : '0'), \
	  (byte & 0x40 ? '1' : '0'), \
	  (byte & 0x20 ? '1' : '0'), \
	  (byte & 0x10 ? '1' : '0'), \
	  (byte & 0x08 ? '1' : '0'), \
	  (byte & 0x04 ? '1' : '0'), \
	  (byte & 0x02 ? '1' : '0'), \
	  (byte & 0x01 ? '1' : '0')

static bool bitstream_pop(struct bitstream *s)
{

	int real_byte = s->byte;
	int real_bit = 7 - s->bit;

	assert(real_bit >= 0);
	assert(real_bit < 8);
	assert(real_byte >= 0);
	assert(real_byte < s->total_bytes);

	if (s->bit == 0)
		printf("\t\tbyte(%d) = %d %x " BYTE_TO_BINARY_PATTERN "\n",
				s->byte,
				s->data[real_byte],
				s->data[real_byte],
				BYTE_TO_BINARY(s->data[real_byte])
		      );

	bool r = s->data[real_byte] & (1 << real_bit);
	s->bit = (s->bit + 1) % 8;
	if (!s->bit)
		s->byte += 1;
	if (s->byte >= s->total_bytes)
		fail("got too much bits from a bitstream!");

	printf("\tbit %d\n", r);
	return r;
}

static unsigned long bitstream_pop_ulong(struct bitstream *s, int nbits)
{
	unsigned long r = 0;
	for (int i = 0; i < nbits; i++)
		r = 2*r + bitstream_pop(s);
	return r;
}

static void bitstream_align_16(struct bitstream *s)
{
	printf("going to align... bit=%d, byte=%d\n", s->bit, s->byte);
	int t[100], cx = 0;
	while (s->bit || s->byte % 2)
		t[cx++] = bitstream_pop(s);
	printf("%d lost bits:", cx);
	for (int i = 0; i < cx; i++)
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

static int global_huffman_tree_brc_0[4][2] = {
	[1] = { 0,  2},
	[2] = {-1,  3},
	[3] = {-2, -3}
};
//       1 --- 2 --- 3 --- -3
//       |     |     |
//       |     |     |
//       0     -1    -2

static int global_huffman_tree_brc_1[5][2] = {
	[1] =  { 0,  2},
	[2] =  {-1,  3},
	[3] =  {-2,  4},
	[4] =  {-3, -4}
};
//       1 --- 2 --- 3 --- 4 --- -4
//       |     |     |     |
//       |     |     |     |
//       0     -1    -2    -3

static int global_huffman_tree_brc_2[7][2] = {
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

static int global_huffman_tree_brc_3[10][2] = {
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
//       |L
//       | L
//       0  -1

static int global_huffman_tree_brc_4[16][2] = {
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

static int (*global_huffman_tree[5])[2] = {
	global_huffman_tree_brc_0,
	global_huffman_tree_brc_1,
	global_huffman_tree_brc_2,
	global_huffman_tree_brc_3,
	global_huffman_tree_brc_4
};

// page 78
static double (*global_table_B_FDBAQ[5]) = {
	[0] = (double[4]){ 3,  3,  3.16,  3.53 },
	[1] = (double[4]){ 4,  4,  4.08,  4.37 },
	[2] = (double[6]){ 6,  6,  6   ,  6.15,  6.5 ,  6.88 },
	[3] = (double[7]){ 9,  9,  9   ,  9   ,  9.36,  9.5 , 10.1},
	[4] = (double[9]){15, 15, 15   , 15   , 15   , 15   ,15.22,15.5,16.05}
};

// page 79
static double (*global_table_NRL_FDBAQ[5]) = {
	[0] = (double[4]){ 0.3637, 1.0915, 1.8208, 2.6406 },
	[1] = (double[5]){ 0.3042, 0.9127, 1.5216, 2.1313, 2.8426 },
	[2] = (double[7]){ 0.2305, 0.6916, 1.1528, 1.6140, 2.0754, 2.5369,
		3.1191 },
	[3] = (double[10]){ 0.1702, 0.5107, 0.8511, 1.1916, 1.5321, 1.8726,
		2.2131, 2.5536, 2.8942, 3.3744 },
	[4] = (double[16]){ 0.1130, 0.3389, 0.5649, 0.7908, 1.0167, 1.2428,
		1.4687, 1.6947, 1.9206, 2.1466, 2.3725, 2.5985, 2.8244, 3.0504,
		3.2764, 3.6623 },
};

// page 80
static double global_table_SF[256] = {
	0.00, 0.63, 1.25, 1.88, 2.51, 3.13, 3.76, 4.39, 5.01, 5.64, 6.27, 6.89,
	7.52, 8.15, 8.77, 9.40, 10.03, 10.65, 11.28, 11.91, 12.53, 13.16,
	13.79, 14.41, 15.04, 15.67, 16.29, 16.92, 17.55, 18.17, 18.80, 19.43,
	20.05, 20.68, 21.31, 21.93, 22.56, 23.19, 23.81, 24.44, 25.07, 25.69,
	26.32, 26.95, 27.57, 28.20, 28.83, 29.45, 30.08, 30.71, 31.33, 31.96,
	32.59, 33.21, 33.84, 34.47, 35.09, 35.72, 36.35, 36.97, 37.60, 38.23,
	38.85, 39.48, 40.11, 40.73, 41.36, 41.99, 42.61, 43.24, 43.87, 44.49,
	45.12, 45.75, 46.37, 47.00, 47.63, 48.25, 48.88, 49.51, 50.13, 50.76,
	51.39, 52.01, 52.64, 53.27, 53.89, 54.52, 55.15, 55.77, 56.40, 57.03,
	57.65, 58.28, 58.91, 59.53, 60.16, 60.79, 61.41, 62.04, 62.98, 64.24,
	65.49, 66.74, 68.00, 69.25, 70.50, 71.76, 73.01, 74.26, 75.52, 76.77,
	78.02, 79.28, 80.53, 81.78, 83.04, 84.29, 85.54, 86.80, 88.05, 89.30,
	90.56, 91.81, 93.06, 94.32, 95.57, 96.82, 98.08, 99.33, 100.58, 101.84,
	103.09, 104.34, 105.60, 106.85, 108.10, 109.35, 110.61, 111.86, 113.11,
	114.37, 115.62, 116.87, 118.13, 119.38, 120.63, 121.89, 123.14, 124.39,
	125.65, 126.90, 128.15, 129.41, 130.66, 131.91, 133.17, 134.42, 135.67,
	136.93, 138.18, 139.43, 140.69, 141.94, 143.19, 144.45, 145.70, 146.95,
	148.21, 149.46, 150.71, 151.97, 153.22, 154.47, 155.73, 156.98, 158.23,
	159.49, 160.74, 161.99, 163.25, 164.50, 165.75, 167.01, 168.26, 169.51,
	170.77, 172.02, 173.27, 174.53, 175.78, 177.03, 178.29, 179.54, 180.79,
	182.05, 183.30, 184.55, 185.81, 187.06, 188.31, 189.57, 190.82, 192.07,
	193.33, 194.58, 195.83, 197.09, 198.34, 199.59, 200.85, 202.10, 203.35,
	204.61, 205.86, 207.11, 208.37, 209.62, 210.87, 212.13, 213.38, 214.63,
	215.89, 217.14, 218.39, 219.65, 220.90, 222.15, 223.41, 224.66, 225.91,
	227.17, 228.42, 229.67, 230.93, 232.18, 233.43, 234.69, 235.94, 237.19,
	238.45, 239.70, 240.95, 242.21, 243.46, 244.71, 245.97, 247.22, 248.47,
	249.73, 250.98, 252.23, 253.49, 254.74, 255.99, 255.99
};


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

// huffman decoder (converts bits into M-codes)
static void extract_scodes(int *scode, struct bitstream *s, int brc, int n)
{
	int (*huf)[2] = global_huffman_tree[brc];
	for (int i = 0; i < n; i++) // flowchart page 70
	{
		int sign = bitstream_pop(s) ? -1 : 1;
		int state = 1;
		while (state > 0)
			state = huf[state][bitstream_pop(s)];
		scode[i] = sign * (-state);
		printf("scode[%d/%d]{%d} = %d\n", i, n, brc, scode[i]);
	}
}

static double sigma_factor(int thidx)
{
	return thidx;
}


static double compute_svalue(int brc, int thidx, int scode)
{
	// page 74
	printf("scode=%d\n", scode);
	int mcode  = abs(scode);
	int sign   = scode<0 ? -1 : 1;
	double B   = global_table_B_FDBAQ[brc][thidx];
	double SF  = global_table_SF[thidx];
	int t[5]   = { 3, 3, 5, 6, 8};   // THIDX threshold for nominal case
	int k[5]   = { 3, 4, 6, 9, 15 }; // mcode threshold inside simple case
	double NRL = global_table_NRL_FDBAQ[brc][mcode];
	if (thidx <= t[brc])
	{
		if (mcode <  k[brc])
			return sign * mcode;
		if (mcode == k[brc])
			return sign * B;
		fail("brc=%d k=%d mcode=%d thidx=%d",brc,k[brc],mcode,thidx);
	}
	return sign * NRL * SF;
}

void s1a_decode_line(complex float *out, struct s1a_isp *x)
{
	s1a_isp_verify_sanity(x);
	s1a_print_some_isp_fields(x);
	s1a_isp_verify_nominal_mode(x);

	// setup init bitstream
	struct bitstream s[1];
	bitstream_init(s, x->data, x->data_size);


	// variables (cf. page 60 of document S1-IF-ASD-PL-0007)
	int NQ = x->secondary_header.field.number_of_quads;
	int NB = ceil(NQ / 128.0); // number of blocks
	//int BRC[NB];               // bit rate of each block
	//int THIDX[NB];             // THIDX of each block
	//int num_hcodes[NB];
	//int code_ie[NQ];
	//int code_io[NQ];
	//int code_qe[NQ];
	//int code_qo[NQ];
	int *BRC        = xmalloc(NB * sizeof(int));
	int *THIDX      = xmalloc(NB * sizeof(int));
	int *num_hcodes = xmalloc(NB * sizeof(int));
	int *code_ie    = xmalloc(NQ * sizeof(int));
	int *code_io    = xmalloc(NQ * sizeof(int));
	int *code_qe    = xmalloc(NQ * sizeof(int));
	int *code_qo    = xmalloc(NQ * sizeof(int));
	for (int i = 0; i < NQ; i++)
		code_ie[i] = code_io[i] = code_qe[i] = code_qo[i] = -42;

	// set number of hcodes for each block
	printf("NB = %d\n", NB);
	for (int b = 0; b < NB; b++)
		num_hcodes[b] = b < NB-1 ? 128 : NQ - 128 * (NB-1); // page 70
	printf("last num_hcodes = %d\n", num_hcodes[NB-1]);

	// decode IE data
	for (int b = 0; b < NB; b++)
	{
		BRC[b] = bitstream_pop_ulong(s, 3);
		printf("scodes IE block b=%d (brc=%d)\n", b, BRC[b]);
		printf("BRC[%d] = %d\n", b, BRC[b]);
		if (BRC[b] < 0 || BRC[b] > 4) fail("bad BRC=%d\n", BRC[b]);

		extract_scodes(code_ie+128*b, s, BRC[b], num_hcodes[b]);
	}
	bitstream_align_16(s);

	// decode IO data
	for (int b = 0; b < NB; b++)
	{
		printf("scodes IO block b=%d\n", b);
		extract_scodes(code_io+128*b, s, BRC[b], num_hcodes[b]);
	}
	bitstream_align_16(s);

	// decode QE data
	for (int b = 0; b < NB; b++)
	{
		THIDX[b] = bitstream_pop_ulong(s, 8);
		printf("scodes QE block b=%d (thidx=%d)\n", b, THIDX[b]);
		printf("THIDX[%d] = %d\n", b, THIDX[b]);

		extract_scodes(code_qe+128*b, s, BRC[b], num_hcodes[b]);
	}
	bitstream_align_16(s);

	// decode QO data
	for (int b = 0; b < NB; b++)
	{
		printf("scodes QO block b=%d\n", b);
		extract_scodes(code_qo+128*b, s, BRC[b], num_hcodes[b]);
	}
	//bitstream_align_16(s);

	// space for quads of svalues
	//double svalue_ie[NQ];
	//double svalue_io[NQ];
	//double svalue_qe[NQ];
	//double svalue_qo[NQ];
	double *svalue_ie = xmalloc(NQ * sizeof(double));
	double *svalue_io = xmalloc(NQ * sizeof(double));
	double *svalue_qe = xmalloc(NQ * sizeof(double));
	double *svalue_qo = xmalloc(NQ * sizeof(double));

	for (int i = 0; i < NQ; i++)
	{
		int b = i/128;
		printf("scode IE (i=%d , b=%d) = %d\n", i, b, code_ie[i]);
	}

	// convert Scodes (small signed integers) to Svalues (physical samples)
	for (int i = 0; i < NQ; i++)
	{
		int b = i / 128; // TODO: verify off-by-one consistency here
		svalue_ie[i] = compute_svalue(BRC[b], THIDX[b], code_ie[i]);
		svalue_io[i] = compute_svalue(BRC[b], THIDX[b], code_io[i]);
		svalue_qe[i] = compute_svalue(BRC[b], THIDX[b], code_qe[i]);
		svalue_qo[i] = compute_svalue(BRC[b], THIDX[b], code_qo[i]);
	}

	for (int i = 0; i < NQ; i++)
	{
		int b = i / 128;
		printf("quad[%d]{%d,%d} = %d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\n", i,
				BRC[b], THIDX[b],
			code_ie[i], code_io[i], code_qe[i], code_qo[i],
			svalue_ie[i], svalue_io[i], svalue_qe[i], svalue_qo[i]
			);
	}

}
