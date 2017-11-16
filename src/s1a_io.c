#include <stdlib.h>
#include <stdio.h>

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "endianness.c"

#include "s1a.h"

// S1-IF-ASD-PL-0007 page 35
#define FILTER_REF_FREQ 37.53472224
static struct { // decimation filter
	int id;
	double decimation_filter_bandwith;
	int decimation_ratio_num;
	int decimation_ratio_den;
	double sampling_frequency_after_decimation;
	int filter_length;
	char *sar_swath;
} s1a_table_decimation_filter[12] = {
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

// S1-IF-ASD-PL-0007 page 19
static const char *s1a_table_ecc_names[48] = {
	[0]  = "contingency (0)", // reserved for ground test or mode upgrading
	[1]  = "stripmap 1",
	[2]  = "stripmap 2",
	[3]  = "stripmap 3",
	[4]  = "stripmap 4",
	[5]  = "stripmap 5-N", // stripmap 5 imaging on northern hemisphere
	[6]  = "stripmap 6",
	[7]  = "contingency (7)", // reserved for ground test or mode upgrading
	[8]  = "interferometric wide swath",
	[9]  = "wave mode", // leapfrog mode
	[10] = "stripmap 5-S", // stripmap 5 imaging on southern hemisphere
	[11] = "stripmap 1 w/0 interl.Cal",
	[12] = "stripmap 2 w/0 interl.Cal",
	[13] = "stripmap 3 w/0 interl.Cal",
	[14] = "stripmap 4 w/0 interl.Cal",
	[15] = "RFC mode", // RF char. mode based on PCC sequences
	[16] = "test mode oper/bypass",
	[17] = "elevation notch S3",
	[18] = "azimuth notch S1",
	[19] = "azimuth notch S2",
	[20] = "azimuth notch S3",
	[21] = "azimuth notch S4",
	[22] = "azimuth notch S5-N",
	[23] = "azimuth notch S5-S",
	[24] = "azimuth notch S6",
	[25] = "stripmap 5-N w/o interl.Cal",
	[26] = "stripmap 5-S w/o interl.Cal",
	[27] = "stripmap 6 w/o interl.Cal",
	[28]  = "contingency (28)", // ground testing or mode upgrading
	[29]  = "contingency (29)", // ground testing or mode upgrading
	[30]  = "contingency (30)", // ground testing or mode upgrading
	[31]  = "elevation notch S3 w/o interl.Cal",
	[32]  = "extra wide swath",
	[33] = "azimuth notch S1 w/o interl.Cal",
	[34] = "azimuth notch S3 w/o interl.Cal",
	[35] = "azimuth notch S6 w/o interl.Cal",
	[36]  = "contingency (36)", // ground testing or mode upgrading
	[37]  = "noise characterisation S1",
	[38]  = "noise characterisation S2",
	[39]  = "noise characterisation S3",
	[40]  = "noise characterisation S4",
	[41]  = "noise characterisation S5-N",
	[42]  = "noise characterisation S5-S",
	[43]  = "noise characterisation S6",
	[44]  = "noise characterisation EWS",
	[45]  = "noise characterisation IWS",
	[46]  = "noise characterisation Wave",
	[47]  = "contingency (47)", // ground testing or mode upgrading
};



static unsigned char xget_byte(FILE *f)
{
	int r = fgetc(f);
	if (r == EOF)
		fail("could not read another byte from file");
	return r;
}

static void xget_bytes(unsigned char *x, FILE *f, int n)
{
	int r = fread(x, 1, n, f);
	if (n != r)
		fail("could not read %d bytes from file (got %d)", n, r);
}

static long get_file_size(FILE *f)
{
	fseek(f, 0, SEEK_END);
	long r = ftell(f);
	rewind(f);
	return r;
}

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
static double s1a_extract_datum_TXPRR(struct s1a_isp *t)
{
	uint16_t x = t->secondary_header.field.tx_ramp_rate;
	int pole   = x >> 15;           // extract first bit
	int txprr  = x & (0xffff >> 1); // remove first bit
	int sign   = pole ? 1 : -1;
	double factor  = FILTER_REF_FREQ * FILTER_REF_FREQ / (1<<21);
	return sign * factor * txprr;
}

// S1-IF-ASD-PL-0007 page 38
static double s1a_extract_datum_TXPSF(struct s1a_isp *t)
{
	uint16_t x = t->secondary_header.field.tx_pulse_start_frequency;
	int pole   = x >> 15;           // extract first bit
	int txpsf  = x & (0xffff >> 1); // remove first bit
	int sign   = pole ? 1 : -1;
	double factor  = FILTER_REF_FREQ * FILTER_REF_FREQ / (1<<14);
	double TXPRR = s1a_extract_datum_TXPRR(t);
	return TXPRR/(4*FILTER_REF_FREQ) + sign * factor * txpsf;
}

// S1-IF-ASD-PL-0007 page 39
static double s1a_extract_datum_TXPL(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.tx_pulse_length;
	return x / FILTER_REF_FREQ;
}

// S1-IF-ASD-PL-0007 page 40
static double s1a_extract_datum_PRI(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.PRI;
	return x / FILTER_REF_FREQ;
}

// S1-IF-ASD-PL-0007 page 41
static double s1a_extract_datum_SWST(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.SWST;
	return x / FILTER_REF_FREQ;
}

// S1-IF-ASD-PL-0007 page 42
static double s1a_extract_datum_SWL(struct s1a_isp *t)
{
	// WARNING: 3 bytes! (TODO : check endianness)
	uint32_t x = t->secondary_header.field.SWL;
	return x / FILTER_REF_FREQ;
}


void s1a_load_whole_datafile_trunc(struct s1a_file *x, char *fname, int max_rec)
{
	FILE *f = xfopen(fname, "r");
	long sf = get_file_size(f);

	// for now, we only read the first ISP
	x->n = 0;
	x->t = xmalloc(10000+(long)(1.1*sf));

	int cx = 0;
	while (1) {
		struct s1a_isp *s = x->t + cx;
		//fprintf(stderr, "ftell(%d) = %zx\n", cx, ftell(f));

		unsigned char header[6];
		for (int i = 0; i < 6; i++)
			header[i] = xget_byte(f);

		s->version_number     = header[0] / 0x40;
		s->id                 = (header[0] % 0x40)*0x100 + header[1];
		s->sequence_control   = header[2]*0x100 + header[3];
		s->packet_data_length = header[4]*0x100 + header[5];
		for (int i = 0; i < 62; i++)
			s->secondary_header.byte[i] = xget_byte(f);

		if (s->packet_data_length > 65534)
			fail("packet len %d too big\n",s->packet_data_length);
		if (s->packet_data_length < 100)
			fail("packet len %d too small\n",s->packet_data_length);

		// table 11, p.61, bottom
		s->data_size = s->packet_data_length + 1 - 62;

		s->data = xmalloc(s->data_size);
		xget_bytes(s->data, f, s->data_size);
		//for (int i = 0; i < s->data_size; i++)
		//	s->data[i] = xget_byte(f);

		cx += 1;
		x->n = cx;

		if (max_rec > 0 && cx > max_rec) break;

		{
			int c = fgetc(f);
			if (c == EOF)
				break;
			else
				ungetc(c, f);
		}
	}

	xfclose(f);

	// correct byte endianness of some fields
	for (int i = 0; i < x->n; i++)
	{
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 31);
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 43);
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 53);
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 54);
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 55);
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 56);
		//reverse_bits_of_byte(x->t[i].secondary_header.byte + 57);
		switch_4endianness(x->t[i].secondary_header.byte +  0, 1);
		switch_2endianness(x->t[i].secondary_header.byte +  4, 1);
		switch_4endianness(x->t[i].secondary_header.byte +  6, 2);
		switch_4endianness(x->t[i].secondary_header.byte + 16, 1);
		switch_2endianness(x->t[i].secondary_header.byte + 21, 1);
		switch_4endianness(x->t[i].secondary_header.byte + 23, 2);
		switch_2endianness(x->t[i].secondary_header.byte + 36, 2);
		switch_3endianness(x->t[i].secondary_header.byte + 40, 1);
		switch_3endianness(x->t[i].secondary_header.byte + 44, 3);
		// TODO: what happens with bytes 54 and 55? it seems ambiguous
		switch_2endianness(x->t[i].secondary_header.byte + 59, 1);
	}

	// TODO: undo sub-commutation of ancillary data (in blocks of 64)
}

void s1a_load_whole_datafile(struct s1a_file *x, char *fname)
{
	s1a_load_whole_datafile_trunc(x, fname, 0);
}

void s1a_load_whole_annot_file(struct s1a_annot_file *x, char *fname)
{
	FILE *f = xfopen(fname, "r");
	int fs = get_file_size(f);
	x->n = fs / 26;
	if (x->n * 26 != fs) fail("bad annotation file size %d\n", fs);
	x->t = xmalloc(x->n * 26);
	int r = fread(x->t, 26, x->n, f);
	if (r != x->n) fail("inconsistent annot file size\n");
	xfclose(f);

	// correct byte endianness of some fields
	for (int i = 0; i < x->n; i++)
	{
		switch_8endianness(x->t[i].record.byte +  0, 2);
		switch_2endianness(x->t[i].record.byte + 16, 3);
	}
}

void s1a_load_whole_index_file(struct s1a_index_file *x, char *fname)
{
	FILE *f = xfopen(fname, "r");
	int fs = get_file_size(f);
	x->n = fs / 36;
	if (x->n * 36 != fs) fail("bad index file size %d\n", fs);
	x->t = xmalloc(x->n * 36);
	int r = fread(x->t, 36, x->n, f);
	if (r != x->n) fail("inconsistent index file size\n");
	xfclose(f);

	// correct byte endianness of some fields
	for (int i = 0; i < x->n; i++)
	{
		switch_8endianness(x->t[i].record.byte +  0, 2);
		switch_4endianness(x->t[i].record.byte + 16, 2);
		switch_8endianness(x->t[i].record.byte + 24, 1);
	}
}

