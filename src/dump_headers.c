#include <stdio.h>

#include "xfopen.c"
#include "xmalloc.c"

#include "s1a.h"

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
	if (i < x->t[j].data_size-1)
		t[j*w+i] = x->t[j].data[i];
		//t[j*w+i] = ((i%2) == 0) ?
		//	hypot(x->t[j].data[i], x->t[j].data[i+1])/1.4
		//	: 0;

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
	//s1a_annot_dump(xa);
	//s1a_index_dump(xi);

	//s1a_dump_image_to_blocks_pgm(filename_out, x);

	return 0;
}
