#include <stdint.h>

// TODO: compile-time asserts to verify exact struct sizes

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
			uint32_t coarse_time                 : 32 ; // 0  _  6
			uint16_t fine_time                   : 16 ; // 4  _ 10

			// fixed ancillary data service (14 bytes, 6-19)
			uint32_t sync_marker                 : 32 ; // 6  _ 12
			uint32_t data_take_id                : 32 ; // 10 _ 16
			uint8_t  ecc_number                  :  8 ; // 14 _ 20
			uint8_t  first_spare_bit             :  1 ; // 15 _ 21
			uint8_t  test_mode                   :  3 ;
			uint8_t  rx_channel_id               :  4 ;
			uint32_t instrument_configuration_id : 32 ; // 16 _ 22

			// sub commutated ancill. data service (3 bytes, 20-22)
			uint8_t  data_word_index             :  8 ; // 20 _ 26
			uint16_t data_word                   : 16 ; // 21 _ 27

			// counter service (8 bytes, 23-30)
			uint32_t space_packet_count          : 32 ; // 23 _ 29
			uint32_t pri_count                   : 32 ; // 27 _ 33

			// radar configuration support service (22 bytes, 31-52)
			//uint8_t  error_flag                  :  1 ; // 31 _ 37
			//uint8_t  zeroth_spare_2bit           :  2 ;
			//uint8_t  baq_mode                    :  5 ;
			//reversed:
			uint8_t  baq_mode                    :  5 ;
			uint8_t  zeroth_spare_2bit           :  2 ;
			uint8_t  error_flag                  :  1 ; // 31 _ 37
			uint8_t  baq_block_length            :  8 ; // 32 _ 38
			uint8_t  spare_byte                  :  8 ; // 33 _ 39
			uint8_t  range_decimation            :  8 ; // 34 _ 40
			uint8_t  rx_gain                     :  8 ; // 35 _ 41
			uint16_t tx_ramp_rate                : 16 ; // 36 _ 42
			uint16_t tx_pulse_start_frequency    : 16 ; // 38 _ 44
			uint32_t tx_pulse_length             : 24 ; // 40 _ 46
			//uint8_t  second_spare_3bit           :  3 ; // 43 _ 49
			//uint8_t  rank                        :  5 ;
			//reversed:
			uint8_t  rank                        :  5 ;
			uint8_t  second_spare_3bit           :  3 ; // 43 _ 49
			uint32_t PRI                         : 24 ; // 44 _ 50
			uint32_t SWST                        : 24 ; // 47 _ 53
			uint32_t SWL                         : 24 ; // 50 _ 56

			// SAS SSB message (3 bytes, 53-55)
			//uint8_t  ssb_flag                    :  1 ; // 53 _ 59
			//uint8_t  polarisation                :  3 ;
			//uint8_t  temperature_compensation    :  2 ;
			//uint8_t  first_spare_2bit            :  2 ;
			//reversed:
			uint8_t  first_spare_2bit            :  2 ;
			uint8_t  temperature_compensation    :  2 ;
			uint8_t  polarisation                :  3 ;
			uint8_t  ssb_flag                    :  1 ; // 53 _ 59
			//uint8_t  elevation_beam_address      :  4 ; // 54 _ 60
			//uint8_t  second_spare_2bit           :  2 ;
			//uint16_t beam_address                : 10 ;
			//reversed:
			uint16_t beam_address                : 10 ;
			uint8_t  second_spare_2bit           :  2 ;
			uint8_t  elevation_beam_address      :  4 ; // 54 _ 60

			// SES SSB message (3 bytes, 56-58)
			//uint8_t  cal_mode                    :  2 ; // 56 _ 60
			//uint8_t  second_spare_bit            :  1 ;
			//uint8_t  tx_pulse_number             :  5 ;
			//reversed:
			uint8_t  tx_pulse_number             :  5 ;
			uint8_t  second_spare_bit            :  1 ;
			uint8_t  cal_mode                    :  2 ; // 56 _ 60
			//uint8_t  signal_type                 :  4 ; // 57 _ 61
			//uint8_t  third_spare_3bit            :  3 ;
			//uint8_t  swap                        :  1 ;
			//reversed:
			uint8_t  swap                        :  1 ;
			uint8_t  third_spare_3bit            :  3 ;
			uint8_t  signal_type                 :  4 ; // 57 _ 61
			uint8_t  swath_number                :  8 ; // 58 _ 62

			// radar sample count service (3 bytes, 59-61)
			uint16_t number_of_quads             : 16 ; // 59 _ 65
			uint8_t  filler_octet                :  8 ; // 61 _ 67
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

// s1a_io.c
void s1a_load_whole_datafile(struct s1a_file *x, char *fname);
void s1a_load_whole_annot_file(struct s1a_annot_file *x, char *fname);
void s1a_load_whole_index_file(struct s1a_index_file *x, char *fname);
void s1a_load_whole_datafile_trunc(struct s1a_file *x, char *fname, int n);

// s1a_decode.c
#include <complex.h>
void s1a_decode_line(complex float *out, struct s1a_isp *x);
