
struct s1a_isp { // Instrument Source Packet

	// packet primary header
	int version_number;     // 3 bits
	int id;                 // 13 bits
	int sequence_control;   // 2 bytes
	int packet_data_length; // 2 bytes

	// packet data field
	unsigned char secondary_header[62];
	unsigned char *data;
	int data_size;


	// the secondary header can be parsed into the following fields
	// by the function "extract_secondary_header"

	// datation service
	unsigned long coarse_time;                 // 4 BYTES
	unsigned long fine_time;                   // 2 BYTES

	// fixed ancillary data service
	unsigned long sync_marker;                 // 4 BYTES
	unsigned long data_take_id;                // 4 BYTES
	unsigned long ecc_number;                  // 1 BYTE
	unsigned long first_spare_bit;             // 1 bit
	unsigned long test_mode;                   // 3 bits
	unsigned long rx_channel_id;               // 4 bits
	unsigned long instrument_configuration_id; // 4 BYTES

	// sub commutation ancillary data service
	unsigned long data_word_index;             // 1 BYTE
	unsigned long data_word;                   // 2 BYTES

	// counter service
	unsigned long space_packet_count;          // 4 BYTES
	unsigned long pri_count;                   // 4 BYTES

	// radar configuration support service
	unsigned long first_spare_3bit;            // 3 bits
	unsigned long baq_mode;                    // 5 bits
	unsigned long baq_block_length;            // 1 BYTE
	unsigned long spare_byte;                  // 1 BYTE
	unsigned long range_decimation;            // 1 BYTE
	unsigned long rx_gain;                     // 1 BYTE
	unsigned long tx_ramp_rate;                // 2 BYTES
	unsigned long tx_pulse_start_frequency;    // 2 BYTES
	unsigned long tx_pulse_length;             // 3 BYTES
	unsigned long second_spare_3bit;           // 3 bits
	unsigned long rank;                        // 5 bits
	unsigned long PRI;                         // 3 BYTES
	unsigned long SWST;                        // 3 BYTES
	unsigned long SWL;                         // 3 BYTES

	// SAS SSB message
	unsigned long ssb_flag;                    // 1 bit
	unsigned long polarisation;                // 3 bits
	unsigned long temperature_compensation;    // 2 bits
	unsigned long first_spare_2bit;            // 2 bits
	unsigned long elevation_beam_address;      // 4 bits
	unsigned long second_spare_2bit;           // 2 bits
	unsigned long beam_address;                // 10 bits

	// SES SSB message
	unsigned long cal_mode;                    // 2 bits
	unsigned long second_spare_bit;            // 1 bit
	unsigned long tx_pulse_number;             // 5 bits
	unsigned long signal_type;                 // 4 bits
	unsigned long third_spare_3bit;            // 3 bits
	unsigned long swap;                        // 1 bit
	unsigned long swath_number;                // 1 BYTE

	// radar sample count service
	unsigned long num_of_quads;                // 2 BYTES
	unsigned long filler_octet;                // 1 BYTE
};

struct s1a_file { // Measurement Data Component, page 64
	int n;
	struct s1a_isp *t;
};

#include <stdlib.h>
#include <stdio.h>

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"

static unsigned char xget_byte(FILE *f)
{
	int r = fgetc(f);
	if (r == EOF)
		fail("could not read another char from file");
	return r;
}

void s1a_load_whole_datafile(struct s1a_file *x, char *fname)
{
	FILE *f = xfopen(fname, "r");

	// for now, we only read the first ISP
	x->n = 1;
	x->t = xmalloc(x->n * sizeof*x->t);

	struct s1a_isp *s = x->t + 0;

	unsigned char header[6];
	for (int i = 0; i < 6; i++)
		header[i] = xget_byte(f);

	s->version_number     = header[0] / 0x40;
	s->id                 = (header[0] % 0x40)*0x100 + header[1];
	s->sequence_control   = header[2]*0x100 + header[3];
	s->packet_data_length = header[4]*0x100 + header[5];
	for (int i = 0; i < 62; i++)
		s->secondary_header[i] = xget_byte(f);

	if (s->packet_data_length > 65534)
		fail("packed data length %d too big\n", s->packet_data_length);
	s->data_size = s->packet_data_length + 1 - 62; // table 11, p.61, bottom

	s->data = xmalloc(s->data_size);
	for (int i = 0; i < s->data_size; i++)
		s->data[i] = xget_byte(f);

	xfclose(f);
}

void s1a_print_info(struct s1a_file *x)
{
	printf("S1A %d ISP\n", x->n);
	for (int i = 0; i < x->n; i++)
	{
		struct s1a_isp *s = x->t + i;
		printf("packet number %d:\n", 1+i);
		printf("\tversion_number     = %d\n", s->version_number);
		printf("\tid                 = %d\n", s->id);
		printf("\tsequence_control   = %d\n", s->sequence_control);
		printf("\tpacket_data_length = %d\n", s->packet_data_length);
	}
}


int main(int c, char *v[])
{
	if (c != 2) return fprintf(stderr, "usage:\n\t%s raw\n", *v);
	char *filename_in = v[1];

	struct s1a_file x[1];
	s1a_load_whole_datafile(x, filename_in);
	s1a_print_info(x);


	return 0;
}
