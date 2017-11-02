#include <stdint.h>

struct s1a_isp { // Instrument Source Packet

	// packet primary header
	int version_number;     // 3 bits
	int id;                 // 13 bits
	int sequence_control;   // 2 bytes
	int packet_data_length; // 2 bytes

	// packet data field
	union { // secondary_header, accessible via bytes or fields
		unsigned char byte[6];
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
		struct {
			uint64_t sensing_time     : 64 ;
			uint64_t downlink_time    : 64 ;
			uint16_t packet_length    : 16 ;
			uint16_t frames           : 16 ;
			uint16_t missing_frames   : 16 ;
			uint8_t  CRC_flag         :  8 ;
			uint8_t  VCID             :  8 ;
			uint8_t  channel          :  8 ;
			uint8_t  spare            :  8 ;

		} field;
	} record;
};

struct s1a_index { // Index Data Component block descriptor, page 67
	union {
		unsigned char byte[36];
		struct {
			uint64_t date_and_time_uint  : 64 ;
			uint64_t delta_time_uint     : 64 ;
			uint32_t delta_size          : 32 ;
			uint32_t data_units_offset   : 32 ;
			uint64_t byte_offset         : 64 ;
			uint8_t  variable_size_flag  :  8 ;
			uint32_t spare               : 24 ;

		} field;
	} record;
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
#define P(x) printf("\t\t" #x " = %u\n", s->secondary_header.field.x)
		printf("\tdatation service:\n");
		P(coarse_time); P(fine_time);
		printf("\tfixed ancillary data service:\n");
		P(sync_marker); P(data_take_id); P(ecc_number);
		P(first_spare_bit);
		P(test_mode); P(rx_channel_id); P(instrument_configuration_id);
		printf("\tsub commutation ancillary data service:\n");
		P(data_word_index); P(data_word);
		printf("\tcounter service:\n");
		P(space_packet_count); P(pri_count);
		printf("\tradar configuration support service:\n");
		P(first_spare_3bit);
		P(baq_mode); P(baq_block_length); P(spare_byte);
		P(range_decimation); P(rx_gain); P(tx_ramp_rate);
		P(tx_pulse_start_frequency); P(tx_pulse_length);
		P(second_spare_3bit);
		P(rank); P(PRI); P(SWST); P(SWL);
		printf("\tSAS SSB message:\n");
		P(ssb_flag); P(polarisation); P(temperature_compensation);
		P(first_spare_2bit);
		P(elevation_beam_address);
		P(second_spare_2bit);
		P(beam_address);
		printf("\tSES SSB message:\n");
		P(cal_mode);
		P(second_spare_bit);
		P(tx_pulse_number); P(signal_type);
		P(third_spare_3bit);
		P(swap); P(swath_number);
		printf("\tradar sample count service:\n");
		P(num_of_quads); P(filler_octet);
	}
}


int main(int c, char *v[])
{
	if (c != 2) return fprintf(stderr, "usage:\n\t%s raw\n", *v);
	char *filename_in = v[1];

	struct s1a_file x[1];
	fprintf(stderr, "sizeof secondary header = %zu\n",
			sizeof(x->t[0].secondary_header));
	s1a_load_whole_datafile(x, filename_in);
	s1a_print_info(x);


	return 0;
}
