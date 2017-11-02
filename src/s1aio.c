
struct s1a_isp { // Instrument Source Packet

	// packet primary header
	int version_number;     // 3 bits
	int id;                 // 13 bits
	int sequence_control;   // 2 bytes
	int packet_data_length; // 2 bytes

	// packet data field
	union { // secondary_header, accessible via bytes or fields
		unsigned char byte[62];
		struct { // field
			// datation service
			unsigned coarse_time                 : 32 ;
			unsigned fine_time                   : 16 ;

			// fixed ancillary data service
			unsigned sync_marker                 : 32 ;
			unsigned data_take_id                : 32 ;
			unsigned ecc_number                  :  8 ;
			unsigned first_spare_bit             :  1 ;
			unsigned test_mode                   :  3 ;
			unsigned rx_channel_id               :  4 ;
			unsigned instrument_configuration_id : 32 ;

			// sub commutation ancillary data service
			unsigned data_word_index             :  8 ;
			unsigned data_word                   : 16 ;

			// counter service
			unsigned space_packet_count          : 32 ;
			unsigned pri_count                   : 32 ;

			// radar configuration support service
			unsigned first_spare_3bit            :  3 ;
			unsigned baq_mode                    :  5 ;
			unsigned baq_block_length            :  8 ;
			unsigned spare_byte                  :  8 ;
			unsigned range_decimation            :  8 ;
			unsigned rx_gain                     :  8 ;
			unsigned tx_ramp_rate                : 16 ;
			unsigned tx_pulse_start_frequency    : 16 ;
			unsigned tx_pulse_length             : 24 ;
			unsigned second_spare_3bit           :  3 ;
			unsigned rank                        :  5 ;
			unsigned PRI                         : 24 ;
			unsigned SWST                        : 24 ;
			unsigned SWL                         : 24 ;

			// SAS SSB message
			unsigned ssb_flag                    :  1 ;
			unsigned polarisation                :  3 ;
			unsigned temperature_compensation    :  2 ;
			unsigned first_spare_2bit            :  2 ;
			unsigned elevation_beam_address      :  4 ;
			unsigned second_spare_2bit           :  2 ;
			unsigned beam_address                : 10 ;

			// SES SSB message
			unsigned cal_mode                    :  2 ;
			unsigned second_spare_bit            :  1 ;
			unsigned tx_pulse_number             :  5 ;
			unsigned signal_type                 :  4 ;
			unsigned third_spare_3bit            :  3 ;
			unsigned swap                        :  1 ;
			unsigned swath_number                :  8 ;

			// radar sample count service
			unsigned num_of_quads                : 16 ;
			unsigned filler_octet                :  8 ;
		} field;
	} secondary_header;
	unsigned char *data;
	int data_size;


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
		s->secondary_header.byte[i] = xget_byte(f);

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
