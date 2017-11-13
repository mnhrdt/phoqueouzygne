#include <stdint.h>

struct s1a_isp { // Instrument Source Packet

	// packet primary header
	int version_number;     // 3 bits
	int id;                 // 13 bits
	int sequence_control;   // 2 bytes
	int packet_data_length; // 2 bytes

	// packet data field
	union { // secondary_header, accessible via bytes or fields
		unsigned char byte[62];
		struct __attribute__((packed)) { // field
			// datation service (6 bytes, 0-5)
			uint32_t coarse_time                 : 32 ; // 0
			uint16_t fine_time                   : 16 ; // 4

			// fixed ancillary data service (14 bytes, 6-19)
			uint32_t sync_marker                 : 32 ; // 6
			uint32_t data_take_id                : 32 ; // 10
			uint8_t  ecc_number                  :  8 ; // 14
			uint8_t  first_spare_bit             :  1 ; // 15
			uint8_t  test_mode                   :  3 ;
			uint8_t  rx_channel_id               :  4 ;
			uint32_t instrument_configuration_id : 32 ; // 16

			// sub commutation ancill. data service (3 bytes, 20-22)
			uint8_t  data_word_index             :  8 ; // 20
			uint16_t data_word                   : 16 ; // 21

			// counter service (8 bytes, 23-30)
			uint32_t space_packet_count          : 32 ; // 23
			uint32_t pri_count                   : 32 ; // 27

			// radar configuration support service (22 bytes, 31-52)
			uint8_t  first_spare_3bit            :  3 ; // 31
			uint8_t  baq_mode                    :  5 ;
			uint8_t  baq_block_length            :  8 ; // 32
			uint8_t  spare_byte                  :  8 ; // 33
			uint8_t  range_decimation            :  8 ; // 34
			uint8_t  rx_gain                     :  8 ; // 35
			uint16_t tx_ramp_rate                : 16 ; // 36
			uint16_t tx_pulse_start_frequency    : 16 ; // 38
			uint32_t tx_pulse_length             : 24 ; // 40
			uint8_t  second_spare_3bit           :  3 ; // 43
			uint8_t  rank                        :  5 ;
			uint32_t PRI                         : 24 ; // 44
			uint32_t SWST                        : 24 ; // 47
			uint32_t SWL                         : 24 ; // 50

			// SAS SSB message (3 bytes, 53-55)
			uint8_t  ssb_flag                    :  1 ; // 53
			uint8_t  polarisation                :  3 ;
			uint8_t  temperature_compensation    :  2 ;
			uint8_t  first_spare_2bit            :  2 ;
			uint8_t  elevation_beam_address      :  4 ; // 54
			uint8_t  second_spare_2bit           :  2 ;
			uint16_t beam_address                : 10 ;

			// SES SSB message (3 bytes, 56-58)
			uint8_t  cal_mode                    :  2 ; // 56
			uint8_t  second_spare_bit            :  1 ;
			uint8_t  tx_pulse_number             :  5 ;
			uint8_t  signal_type                 :  4 ; // 57
			uint8_t  third_spare_3bit            :  3 ;
			uint8_t  swap                        :  1 ;
			uint8_t  swath_number                :  8 ; // 58

			// radar sample count service (3 bytes, 59-61)
			uint16_t num_of_quads                : 16 ; // 59
			uint8_t  filler_octet                :  8 ; // 61
		} field;
	} secondary_header;
	unsigned char *data;
	int data_size;
};

struct s1a_file { // Measurement Data Component, page 64
	int n;
	struct s1a_isp *t;
};

struct s1a_annot { // Annotation Data Component record, page 65
	union {
		unsigned char byte[26];
		struct __attribute__((packed)) {
			uint64_t sensing_time     : 64 ; // 0
			uint64_t downlink_time    : 64 ; // 8
			uint16_t packet_length    : 16 ; // 16
			uint16_t frames           : 16 ; // 18
			uint16_t missing_frames   : 16 ; // 20
			uint8_t  CRC_flag         :  8 ; // 22
			uint8_t  VCID             :  8 ; // 23
			uint8_t  channel          :  8 ; // 24
			uint8_t  spare            :  8 ; // 25

		} field;
	} record;
};

struct s1a_index { // Index Data Component block descriptor, page 67
	union {
		unsigned char byte[36];
		struct __attribute__((packed)) {
			//uint64_t date_and_time  : 64 ; // 0
			//uint64_t delta_time     : 64 ; // 8
			double   date_and_time;             // 0
			double   delta_time;                // 8
			uint32_t delta_size          : 32 ; // 16
			uint32_t data_units_offset   : 32 ; // 20
			uint64_t byte_offset         : 64 ; // 24
			uint8_t  variable_size_flag  :  8 ; // 32
			uint32_t spare               : 24 ; // 33

		} field;
	} record;
};

struct s1a_annot_file {
	int n;
	struct s1a_annot *t;
};

struct s1a_index_file {
	int n;
	struct s1a_index *t;
};

#include <stdlib.h>
#include <stdio.h>

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "endianness.c"

static unsigned char xget_byte(FILE *f)
{
	int r = fgetc(f);
	if (r == EOF)
		fail("could not read another char from file");
	return r;
}

static long get_file_size(FILE *f)
{
	fseek(f, 0, SEEK_END);
	long r = ftell(f);
	rewind(f);
	return r;
}


void s1a_load_whole_datafile(struct s1a_file *x, char *fname)
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
		for (int i = 0; i < s->data_size; i++)
			s->data[i] = xget_byte(f);

		cx += 1;
		x->n = cx;

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

void s1a_dump_headers(struct s1a_file *x)
{
#define P(f) printf("mdc[%d]." #f " = %lu\n", i, (unsigned long)x->t[i].f)
#define P2(f) printf("mdc[%d]." #f " = %lu\n", i, \
		(unsigned long)x->t[i].secondary_header.field.f)
	for (int i = 0; i < x->n; i++)
	{
		P(version_number);
		P(id);
		P(sequence_control);
		P(packet_data_length);

		P2(coarse_time); P2(fine_time); P2(sync_marker);
		P2(data_take_id); P2(ecc_number); P2(first_spare_bit);
		P2(test_mode); P2(rx_channel_id);
		P2(instrument_configuration_id); P2(data_word_index);
		P2(data_word); P2(space_packet_count); P2(pri_count);
		P2(first_spare_3bit); P2(baq_mode); P2(baq_block_length);
		P2(spare_byte); P2(range_decimation); P2(rx_gain);
		P2(tx_ramp_rate); P2(tx_pulse_start_frequency);
		P2(tx_pulse_length); P2(second_spare_3bit); P2(rank); P2(PRI);
		P2(SWST); P2(SWL); P2(ssb_flag); P2(polarisation);
		P2(temperature_compensation); P2(first_spare_2bit);
		P2(elevation_beam_address); P2(second_spare_2bit);
		P2(beam_address); P2(cal_mode); P2(second_spare_bit);
		P2(tx_pulse_number); P2(signal_type); P2(third_spare_3bit);
		P2(swap); P2(swath_number); P2(num_of_quads); P2(filler_octet);
	}
#undef P
#undef P2
}

void s1a_annot_dump(struct s1a_annot_file *x)
{
#define P(f) printf("annot[%d]." #f " = %lu\n", \
	       i, (unsigned long)x->t[i].record.field.f)
	for (int i = 0; i < x->n; i++)
	{
		P(sensing_time);
		P(downlink_time);
		P(packet_length);
		P(frames);
		P(missing_frames);
		P(CRC_flag);
		P(VCID);
		P(channel);
		P(spare);
	}
#undef P
}

static void pgm_write(char *fname, uint8_t *x, int w, int h)
{
	FILE *f = xfopen(fname, "w");
	fprintf(f, "P5\n%d %d\n255\n", w, h);
	fwrite(x, h, w, f);
	xfclose(f);
}

void s1a_dump_image_to_blocks_pgm(char *fname, struct s1a_file *x)
{
	int h = x->n;
	int w = 0;
	for (int i = 0; i < h; i++)
		if (x->t[i].data_size > w)
			w = x->t[i].data_size;
	uint8_t *t = xmalloc(w*h);
	memset(t, 0, w*h);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	if (i < x->t[j].data_size)
		t[j*w+i] = x->t[j].data[i];

	pgm_write(fname, t, w, h);
}


void s1a_index_dump(struct s1a_index_file *x)
{
#define P(f) printf("index[%d]." #f " = %lu\n", \
	       i, (unsigned long)x->t[i].record.field.f)
	for (int i = 0; i < x->n; i++)
	{
		P(date_and_time);
		P(delta_time);
		P(delta_size);
		P(data_units_offset);
		P(byte_offset);
		P(variable_size_flag);
		P(spare);
	}
#undef P
}

// (from page 119)
//
// The Measurement Data Component is a single binary file, including X-Band
// Recorded HKTM data in the form of transfer frames according to the following
// specification.
//
// * The binary file includes HKTM transfer frames as received in X-Band from
//   the 4 on-board memory Packet Stores allocated to HKTM storage.
//
// * The transfer frames sequence is as per Packet Stores downlink order,
//   starting with the PS with associated lower VC ID.
//
// * Transfer frames do not include the Attached Synchronisation Marker (ASM)
//   and the Reed-Solomon (R-S) code block (a single transfer frame is 1912
//   octets).
//
// * Transfer frames are concatenated with no separator between consecutive
//   frames.
//
// * Transfer frames are not randomised.
//
// * Only transfer frames that have passed the R-S decoding are included in the
//   binary file.

#include <assert.h>

int main(int c, char *v[])
{
	if (c != 5) return fprintf(stderr, "usage:\n\t"
	             "%s in-raw in-annot in-index out.pgm >meta.txt\n", *v);
	//             0 1      2        3        4
	char *filename_x   = v[1];
	char *filename_xa  = v[2];
	char *filename_xi  = v[3];
	char *filename_out = v[4];

	struct s1a_file        x[1];
	struct s1a_annot_file xa[1];
	struct s1a_index_file xi[1];

	assert(62 == sizeof x->t[0].secondary_header);
	assert(26 == sizeof xa->t[0]);
	assert(36 == sizeof xi->t[0]);

	s1a_load_whole_datafile  (x , filename_x );
	s1a_load_whole_annot_file(xa, filename_xa);
	s1a_load_whole_index_file(xi, filename_xi);

	printf("%s: %d records\n", filename_x , x ->n);
	printf("%s: %d records\n", filename_xa, xa->n);
	printf("%s: %d records\n", filename_xi, xi->n);

	s1a_dump_headers(x);
	s1a_annot_dump(xa);
	s1a_index_dump(xi);

	s1a_dump_image_to_blocks_pgm(filename_out, x);

	return 0;
}
