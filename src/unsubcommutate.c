#include <stdio.h>

#include "xfopen.c"
#include "xmalloc.c"
#include "bits.c"

#include "s1a.h"

static void dump_8_bits(FILE *f, uint8_t x)
{
	fprintf(stderr, "%d:" BYTE_TO_BINARY_PATTERN, x, BYTE_TO_BINARY(x));
}
static void dump_16_bits(FILE *f, uint16_t x)
{
	fprintf(stderr, "%d:"BYTE_TO_BINARY_PATTERN","BYTE_TO_BINARY_PATTERN,
			x, BYTE_TO_BINARY(x>>8),BYTE_TO_BINARY(x));
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
		P2(error_flag);
		P2(zeroth_spare_2bit);
		P2(baq_mode); P2(baq_block_length);
		P2(spare_byte); P2(range_decimation); P2(rx_gain);
		P2(tx_ramp_rate); P2(tx_pulse_start_frequency);
		P2(tx_pulse_length); P2(second_spare_3bit); P2(rank); P2(PRI);
		P2(SWST); P2(SWL); P2(ssb_flag); P2(polarisation);
		P2(temperature_compensation); P2(first_spare_2bit);
		P2(elevation_beam_address); P2(second_spare_2bit);
		P2(beam_address); P2(cal_mode); P2(second_spare_bit);
		P2(tx_pulse_number); P2(signal_type); P2(third_spare_3bit);
		P2(swap); P2(swath_number); P2(number_of_quads);
		P2(filler_octet);

		if (1 == x->t[i].secondary_header.field.data_word_index)
		{
			struct s1a_unsubcommutated_block *S = x->t[i].SC;
			printf("SC[%d].dummy = %d\n", i, S->block.field.dummy);
			printf("SC[%d].position = %g %g %g\n", i,
					S->block.field.position[0],
					S->block.field.position[1],
					S->block.field.position[2]);
			printf("SC[%d].velocity = %g %g %g\n", i,
					S->block.field.velocity[0],
					S->block.field.velocity[1],
					S->block.field.velocity[2]);
			printf("SC[%d].GPS_time_POD = %d %d %d %d\n", i,
					S->block.field.GPS_time_POD[0],
					S->block.field.GPS_time_POD[1],
					S->block.field.GPS_time_POD[2],
					S->block.field.GPS_time_POD[3]);
			printf("SC[%d].Q = %g %g %g %g\n", i,
					S->block.field.Q[0],
					S->block.field.Q[1],
					S->block.field.Q[2],
					S->block.field.Q[3]);
			printf("SC[%d].w = %g %g %g\n", i,
					S->block.field.w[0],
					S->block.field.w[1],
					S->block.field.Q[2]);
			printf("SC[%d].GPS_time_DTS = %d %d %d %d\n", i,
					S->block.field.GPS_time_DTS[0],
					S->block.field.GPS_time_DTS[1],
					S->block.field.GPS_time_DTS[2],
					S->block.field.GPS_time_DTS[3]);
			//printf("SC[%d].pointing_status = %d\n", i,
			//		S->block.field.pointing_status);
			printf("SC[%d].pointing_status2 = %d:"BYTE_TO_BINARY_PATTERN","BYTE_TO_BINARY_PATTERN"\n",
			i, S->block.field.pointing_status,
			BYTE_TO_BINARY(S->block.field.pointing_status >> 8),
			BYTE_TO_BINARY(S->block.field.pointing_status)
					);
			printf("SC[%d].temp_status = %d\n", i,
					S->block.field.temp_status);
			// TODO (?): print all temp tiles
			printf("SC[%d].temp_TGU = %d\n", i,
					S->block.field.temp_TGU);

		}
	}
#undef P
#undef P2
}

// fill-in the subcommutated ancillary data
static void s1a_unsubcommutate(struct s1a_file *x)
{
	for (int i = 0; i < x->n; i++)
	{
		int n = -1; // number of valid words so far

		struct s1a_isp *t = x->t + i;
		if (1 == t->secondary_header.field.data_word_index)
		{
			n = 0;

			// first fill-in the whole data of this line
			for (int j = 0; j < 64; j++)
			{
				struct s1a_isp *T = x->t + i + j;
				// extract word and index
				int w = T->secondary_header.field.data_word;
				int y = T->secondary_header.field.data_word_index;

				// check consistency, otherwise restart cycle
				if (y != j+1)
					break;

				// put data on unsubcommutated block
				t->SC->block.byte[2*y+0] = w / 256;
				t->SC->block.byte[2*y+1] = w % 256;
				n += 1;
			}
			// Table 3.2-5 page 23
			switch_8endianness(t->SC->block.byte +  2, 3);
			switch_4endianness(t->SC->block.byte + 26, 3);
			// Table 3.2-6 page 24
			switch_4endianness(t->SC->block.byte + 46, 7);

			// then propagate to the next 64 lines
			fprintf(stdout, "n = %d\n", n);

		}
		//else fprintf(stderr, "%d ignored\n", t->secondary_header.field.data_word_index);
		// TODO: study whether using "previous" or "centered" makes any
		// meaningful difference here
	}
}


#include <assert.h>

int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr,"usage:\n\t%s raw.dat >fields.txt\n",*v);
	//                                        0 1
	char *filename_x = v[1];

	struct s1a_file x[1];
	assert(62 == sizeof x->t[0].secondary_header);

	{
		struct s1a_unsubcommutated_block s[1];
		printf("sizeof s = %zu\n", sizeof *s);
		printf("sizeof s.block = %zu\n", sizeof s->block);
		printf("sizeof s.block = %zu\n", sizeof s->block);
		printf("sizeof s.block.field = %zu\n", sizeof s->block.field);
	}



	s1a_load_whole_datafile  (x , filename_x );

	printf("%s: %d records\n", filename_x , x ->n);

	s1a_unsubcommutate(x);

	s1a_dump_headers(x);
	return 0;
}
