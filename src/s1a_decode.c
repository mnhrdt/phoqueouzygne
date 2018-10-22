#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "fail.c"
#include "bits.c"

#include "s1a.h"

// data structure for bitwise memory traversal
struct bitstream {
	int byte;
	int bit;
	int total_bytes;
	uint8_t *data;
};

// setup a bitstream from a memory array
static void bitstream_init(struct bitstream *s, void *data, int total_bytes)
{
	s->byte = 0;
	s->bit = 0;
	s->data = data;
	s->total_bytes = total_bytes;
}

// extract next bit from the bitstream
static bool bitstream_pop(struct bitstream *s)
{
	// check if there is still available data
	if (s->byte >= s->total_bytes)
		fail("got too much bits from a bitstream!");

	// extract bit
	bool r = s->data[s->byte] & (1 << (7 - s->bit));

	// increase counters
	s->bit = (s->bit + 1) % 8;
	if (!s->bit)
		s->byte += 1;

	return r;
}

// extract nbits to form an integer
static unsigned long bitstream_pop_ulong(struct bitstream *s, int nbits)
{
	unsigned long r = 0;
	for (int i = 0; i < nbits; i++)
		r = 2*r + bitstream_pop(s);
	return r;
}

// extract enough bits so that the head is aligned to a multiple of 16
static void bitstream_align_16(struct bitstream *s)
{
	int t[100], cx = 0;
	while (s->bit || s->byte % 2)
		t[cx++] = bitstream_pop(s);
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



// verify that some constant fields have the appropriate constants
void s1a_isp_verify_sanity(struct s1a_isp *x)
{
	if (x->secondary_header.field.sync_marker != 0x352ef853)
		fail("incorrect sync marker found");

	if (x->secondary_header.field.baq_block_length != 31)
		fail("incorrect baq_block_field");
}

// check wether the line is a nominal datum (FDBAQ-encoded echo)
bool s1a_isp_nominal_mode_P(struct s1a_isp *x)
{
	int TSTMOD = x->secondary_header.field.test_mode;
	int SIGTYP = x->secondary_header.field.signal_type;
	int BAQMOD = x->secondary_header.field.baq_mode;
	return   TSTMOD == 0   &&   SIGTYP == 0   &&   BAQMOD == 12;
	//if (TSTMOD !=  0) fail("bad TSTMOD %d (should be 0)\n", TSTMOD);
	//if (SIGTYP !=  0) fail("bad SIGTYP %d (should be 0)\n", SIGTYP);
	//if (BAQMOD != 12) fail("bad BAQMOD %d (should be 12)\n", BAQMOD);
}

// debugging function to dump most data fields of the header
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

// huffman decoder (converts bits into S-codes)
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
	}
}

// compute a physical S-value (double) from a quantized S-code (int)
// S1-IF-ASD-PL-0007 page 74
static double compute_svalue(int brc, int thidx, int scode)
{
	// page 74
	int mcode  = abs(scode);
	int sign   = scode<0 ? -1 : 1;
	double B   = global_table_B_FDBAQ[brc][thidx];
	double SF  = global_table_SF[thidx];
	double NRL = global_table_NRL_FDBAQ[brc][mcode];
	int t[5]   = { 3, 3, 5, 6, 8};   // THIDX threshold for nominal case
	int k[5]   = { 3, 4, 6, 9, 15 }; // mcode threshold inside simple case
	if (thidx <= t[brc])
	{
		if (mcode <  k[brc])
			return sign * mcode;
		if (mcode == k[brc])
			return sign * B;
		fail("brc=%d k=%d mcode=%d thidx=%d", brc,k[brc],mcode,thidx);
	}
	return sign * NRL * SF;
}

// decode a line into an array of complex numbers (cf. flowchart on page 69)
// returns the number of complex samples filled-in
int s1a_decode_line_fancy(
		complex float *out,
		uint8_t *out_block,
		uint8_t *out_brc,
		uint8_t *out_thidx,
		struct s1a_isp *x)
{
	s1a_isp_verify_sanity(x);
	if (!s1a_isp_nominal_mode_P(x))
		return 0 * fprintf(stderr, "isp %p is non-nominal\n", (void*)x);

	// setup init bitstream
	struct bitstream s[1];
	bitstream_init(s, x->data, x->data_size);

	// variables (cf. page 60 of document S1-IF-ASD-PL-0007)
	int NQ = x->secondary_header.field.number_of_quads;
	int NB = ceil(NQ / 128.0); // number of blocks
	int BRC[NB];               // bit rate of each block
	int THIDX[NB];             // THIDX of each block
	int num_hcodes[NB];
	int code_ie[NQ], code_io[NQ], code_qe[NQ], code_qo[NQ];

	for (int b = 0; b < NB; b++) // set number of hcodes for each block
		num_hcodes[b] = b < NB-1 ? 128 : NQ - 128 * (NB-1); // page 70

	for (int b = 0; b < NB; b++) // decode IE data
	{
		BRC[b] = bitstream_pop_ulong(s, 3);
		if (BRC[b] < 0 || BRC[b] > 4) fail("bad BRC=%d\n", BRC[b]);
		//fprintf(stderr, "BRC[%d] = %d\n", b, BRC[b]);
		extract_scodes(code_ie+128*b, s, BRC[b], num_hcodes[b]);
	}
	bitstream_align_16(s);

	for (int b = 0; b < NB; b++) // decode IO data
		extract_scodes(code_io+128*b, s, BRC[b], num_hcodes[b]);
	bitstream_align_16(s);

	for (int b = 0; b < NB; b++) // decode QE data
	{
		THIDX[b] = bitstream_pop_ulong(s, 8);
		//fprintf(stderr, "THIDX[%d] = %d\n", b, THIDX[b]);
		extract_scodes(code_qe+128*b, s, BRC[b], num_hcodes[b]);
	}
	bitstream_align_16(s);

	for (int b = 0; b < NB; b++) // decode QO data
		extract_scodes(code_qo+128*b, s, BRC[b], num_hcodes[b]);
	bitstream_align_16(s);

	// space for quads of svalues
	double svalue_ie[NQ], svalue_io[NQ], svalue_qe[NQ], svalue_qo[NQ];

	for (int i = 0; i < NQ; i++) // convert Scodes to Svalues
	{
		int b = i / 128; // TODO: verify off-by-one consistency here
		svalue_ie[i] = compute_svalue(BRC[b], THIDX[b], code_ie[i]);
		svalue_io[i] = compute_svalue(BRC[b], THIDX[b], code_io[i]);
		svalue_qe[i] = compute_svalue(BRC[b], THIDX[b], code_qe[i]);
		svalue_qo[i] = compute_svalue(BRC[b], THIDX[b], code_qo[i]);
	}

	for (int i = 0; i < NQ; i++) // interlace even and odd samples
	{
		out[2*i+0] = svalue_ie[i] + I * svalue_qe[i];
		out[2*i+1] = svalue_io[i] + I * svalue_qo[i];
	}

	for (int i = 0; i < NQ; i++)
	{
		int b = i / 128;
		if (out_block) out_block[2*i] = out_block[2*i+1] = b;
		if (out_brc  ) out_brc  [2*i] = out_brc  [2*i+1] = BRC[b];
		if (out_thidx) out_thidx[2*i] = out_thidx[2*i+1] = THIDX[b];
	}

	return 2*NQ;
}

// decode a line into an array of complex numbers (cf. flowchart on page 69)
// returns the number of complex samples filled-in
int s1a_decode_line(complex float *out, struct s1a_isp *x)
{
	return s1a_decode_line_fancy(out, NULL, NULL, NULL, x);
}



///////////////////////////
///////////////////////////
///////////////////////////


// utility functions to extract data with the correct units/conversions

// S1-IF-ASD-PL-0007 page 35
static struct { // decimation filter
	int id;
	double decimation_filter_bandwith;
	int decimation_ratio_num;
	int decimation_ratio_den;
	double sampling_frequency_after_decimation;
	int filter_length;
	char *sar_swath;
} global_table_decimation_filter[12] = {
	[0]  = {0, 100.00, 3,4,  3*4*FILTER_REF_FREQ/4 , 28, "full bandwith" },
	[1]  = {1,  87.71, 2,3,  2*4*FILTER_REF_FREQ/3 , 28, "s1,wv1" },
	[2]  = {2,  0    , 0,0,  0,                       0,   "" },
	[3]  = {3,  74.25, 5,9,  5*4*FILTER_REF_FREQ/9 , 32, "s2" },
	[4]  = {4,  59.44, 4,9,  4*4*FILTER_REF_FREQ/9 , 40, "s3" },
	[5]  = {5,  50.62, 3,8,  3*4*FILTER_REF_FREQ/8 , 48, "s4" },
	[6]  = {6,  44.89, 1,3,  1*4*FILTER_REF_FREQ/3 , 52, "s5" },
	[7]  = {7,  22.2 , 1,6,  1*4*FILTER_REF_FREQ/6 , 92, "ew1" },
	[8]  = {8,  56.59, 3,7,  3*4*FILTER_REF_FREQ/7 , 36, "iw1" },
	[9]  = {9,  42.86, 5,16, 5*4*FILTER_REF_FREQ/16, 68, "s6,iw3" },
	[10] = {10, 15.1 , 3,26, 3*4*FILTER_REF_FREQ/26,120, "ew2,ew3,ew4,ew5"},
	[11] = {11, 48.35, 4,11, 4*4*FILTER_REF_FREQ/11, 44, "iw1,wv2" }
};

// S1-IF-ASD-PL-007 page 77 table 5.1-2
static int global_table_FOO[17] = {
	[0] = 87,
	[1] = 87,
	[2] = 0,
	[3] = 88,
	[4] = 90,
	[5] = 92,
	[6] = 93,
	[7] = 103,
	[8] = 89,
	[9] = 97,
	[10] = 110,
	[11] = 91,
	[12] = 0,
	[13] = 0,
	[14] = 0,
	[15] = 0,
	[16] = 0
};

// S1-IF-ASD-PL-007 page 76 table 5.1-1
#define x -1
static double global_table_CD[12][26] = {
	//     0 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 20 1 2 3 4 5
	[0] = {1,1,2,3,x,x,x,x,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[1] = {1,1,2,x,x,x,x,x,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[2] = {x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[3] = {1,1,2,2,3,3,4,4,5,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[4] = {0,1,1,2,2,3,3,4,4,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[5] = {0,1,1,1,2,2,3,3,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[6] = {0,0,1,x,x,x,x,x,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[7] = {0,0,0,0,0,1,x,x,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[8] = {0,1,1,2,2,3,3,x,x,x, x,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x},
	[9] = {0,0,1,1,1,2,2,2,2,3, 3,3,4,4,4,5,x,x,x,x, x,x,x,x,x,x},
	[10]= {0,0,0,0,0,0,0,1,1,1, 1,1,1,1,1,1,2,2,2,2, 2,2,2,2,3,3},
	[11]= {0,1,1,1,2,2,3,3,3,4, 4,x,x,x,x,x,x,x,x,x, x,x,x,x,x,x}
};
#undef x

static double s1a_extract_datum_TFINE(struct s1a_isp *t)
{
	return (t->secondary_header.field.fine_time + 0.5) / 0x10000;
}

// S1-IF-ASD-PL-0007 page 36
static double s1a_extract_datum_RXG(struct s1a_isp *t)
{
	return -0.5 * t->secondary_header.field.rx_gain;
}

// S1-IF-ASD-PL-0007 page 37
double s1a_extract_datum_TXPRR(struct s1a_isp *t)
{
	uint16_t x = t->secondary_header.field.tx_ramp_rate;
	int pole   = x >> 15;           // extract first bit
	int txprr  = x & (0xffff >> 1); // remove first bit
	int sign   = pole ? 1 : -1;
	double factor  = FILTER_REF_FREQ * FILTER_REF_FREQ / (1<<21);
	return sign * factor * txprr;// * 1e12;
}

// S1-IF-ASD-PL-0007 page 38
double s1a_extract_datum_TXPSF(struct s1a_isp *t)
{
	uint16_t x = t->secondary_header.field.tx_pulse_start_frequency;
	int pole   = x >> 15;           // extract first bit
	int txpsf  = x & (0xffff >> 1); // remove first bit
	int sign   = pole ? 1 : -1;
	double factor  = FILTER_REF_FREQ / (1<<14);
	double TXPRR = s1a_extract_datum_TXPRR(t);
	return /*1e6*(1e-12*/(TXPRR/(4*FILTER_REF_FREQ) + sign * factor * txpsf);
}

// S1-IF-ASD-PL-0007 page 39
double s1a_extract_datum_TXPL(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.tx_pulse_length;
	return (x / FILTER_REF_FREQ);// * 1e-6;
}

// width of the filter (chirp)
int s1a_extract_datum_TXPL3(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.tx_pulse_length;
	if (x < 128 || x > 4223)
		fprintf(stderr, "WARNING: TXPL = %d\n", x);
	int i = t->secondary_header.field.range_decimation;;
	if (i < 0 || i > 11)
		fprintf(stderr, "WARNING: RGDEC = %d\n", i);
	double fdec = global_table_decimation_filter[i].sampling_frequency_after_decimation;
	return ceil(fdec * x / FILTER_REF_FREQ);
}

// width of the filter (chirp)
int s1a_extract_datum_NF(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.tx_pulse_length;
	if (x < 128 || x > 4223)
		fprintf(stderr, "WARNING: TXPL = %d\n", x);
	int i = t->secondary_header.field.range_decimation;;
	if (i < 0 || i > 11)
		fprintf(stderr, "WARNING: RGDEC = %d\n", i);
	return global_table_decimation_filter[i].filter_length;
}

int s1a_extract_datum_TXPL1(struct s1a_isp *t)
{
	return 8 * t->secondary_header.field.tx_pulse_length;
}

int s1a_extract_datum_TXPL2(struct s1a_isp *t)
{
	return 4 * t->secondary_header.field.tx_pulse_length;
}

// S1-IF-ASD-PL-0007 page 40
double s1a_extract_datum_PRI(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.PRI;
	return x / FILTER_REF_FREQ;
}

// S1-IF-ASD-PL-0007 page 41
double s1a_extract_datum_SWST(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.SWST;
	return x / FILTER_REF_FREQ;
}

// S1-IF-ASD-PL-0007 page 42
double s1a_extract_datum_SWL(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.SWL;
	return x / FILTER_REF_FREQ;
}

// S1-IF-ASD-PL-0007 page 42, first row of table
double s1a_extract_datum_SWL1(struct s1a_isp *t)
{
	double x = t->secondary_header.field.SWL;
	return 8 * x;
}

// S1-IF-ASD-PL-0007 page 42, second row of table
double s1a_extract_datum_SWL2(struct s1a_isp *t)
{
	double x = t->secondary_header.field.SWL;
	return 4 * x;
}

// S1-IF-ASD-PL-0007 page 42, third row of table
int s1a_extract_datum_SWL3(struct s1a_isp *t)
{
	int i = t->secondary_header.field.range_decimation; // filter number
	assert(i >= 0);
	assert(i <= 11);
	assert(i != 2);
	int SWL = t->secondary_header.field.SWL;
	int FOO = global_table_FOO[i];
	assert(FOO > 0);
	int B = 2 * SWL - FOO - 17;
	int L = global_table_decimation_filter[i].decimation_ratio_num;
	int M = global_table_decimation_filter[i].decimation_ratio_den;
	int C = B - M * (B/M);
	assert(C >= 0);
	assert(C < 26);
	int D = global_table_CD[i][C];
	assert(D >= 0);
	assert(D <= 5);
	return 2*(L*(B/M) + D + 1);
}



