/*
 ============================================================================
 Name        : compress.c
 Author      : nanash1
 Version     : 0.95
 Description : Compresses Grandia Saturn graphics
 ============================================================================
 */

/**
* This software is supplied "AS IS" without any warranties of
* any kind. The author, nanash1, assumes no responsibility
* or liability for the use of the software, conveys no license or rights under any
* patent, copyright, mask work right, or any other intellectual property rights in
* or to any products. The author also makes no
* representation or warranty that such application will be suitable for the
* specified use without further testing or modification.
*
* Permission to use, copy, modify, and distribute this software and its
* documentation is hereby granted. This copyright, permission, and disclaimer
* notice must appear in all copies of this code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <bits/byteswap.h>
#include "compression.h"

static uint8_t compr_lut_shifts[16];
static uint8_t compr_lut_offset[16];
static uint8_t compr_shifts_stop;
static uint8_t compr_offset_stop;

/**
 * @brief	Generates the decompression lut
 * @param	p_comp_lut		Pointer to the memory the table is written at
 * @return	Status 0:Success 1:Error
 */
static uint8_t generate_lut(uint8_t* p_comp_lut)
{
	uint8_t rtn = 0;
	uint8_t bit_num;
	uint8_t max_bit_num;
	uint8_t decode_char;
	uint16_t offset = 0;

	for(uint_fast8_t j = 0; j < 17; j++){
		max_bit_num = 0;
		for(uint_fast8_t i = 0; i < 17; i++){
			if (i%2){
				bit_num = *(p_comp_lut + (i/2)) & 0xf;
			} else {
				bit_num = *(p_comp_lut + (i/2)) >> 4;
			}

			if (bit_num > max_bit_num){
				max_bit_num = bit_num;
				decode_char = i;
			}
		}

		if (max_bit_num > 8){
			rtn = 1;
			break;
		} else if (max_bit_num == 0){
			rtn = 0;
			break;
		}

		if (decode_char%2){
			*(p_comp_lut + (decode_char/2)) &= 0xf0;
		} else {
			*(p_comp_lut + (decode_char/2)) &= 0xf;
		}

		if (decode_char == 16){
			compr_offset_stop = offset;
			compr_shifts_stop = max_bit_num;
		} else {
			compr_lut_offset[decode_char] = offset;
			compr_lut_shifts[decode_char] = max_bit_num;
		}
		offset += lut_seed[max_bit_num];
	}

	/* check if the table is within bounds */
	if (offset != 256){
		rtn = 1;
	}
	return rtn;

}

/* input data variables */
uint8_t* p_input_data;
uint_fast8_t input_bit_pos = 0;

/* RLE output variables */
uint8_t* p_rle_data;
uint_fast8_t rle_bit_pos = 1;

/* table output variables */
uint8_t* p_delta_data;
uint_fast8_t delta_bit_pos = 0;

uint8_t* p_end_of_file;
uint8_t end_of_file = 0;

/**
 * @brief	Stores or restores the current position in the input stream
 * @param	restore_flag		If set, position is restored
 * @return	Nothing
 */
static inline void store_input_pos(uint8_t restore_flag)
{
	static uint8_t* last_p_input_data;
	static uint_fast8_t last_input_bit_pos;

	if (restore_flag){
		p_input_data = last_p_input_data;
		input_bit_pos = last_input_bit_pos;
	} else {
		last_p_input_data = p_input_data;
		last_input_bit_pos = input_bit_pos;
	}
}

/**
 * @brief	Retrieves the next 4 bits of data from the input stream
 * @return	Next 4 bits of input stream data
 */
static uint8_t get_next_input_bits(void)
{
	uint8_t rtn;

	if (p_input_data > p_end_of_file){
		end_of_file = 1;
		return 0;
	}

	rtn = (*p_input_data >> (4 - input_bit_pos)) & 0xf;
	input_bit_pos += 4;
	if (input_bit_pos > 4){
		input_bit_pos = 0;
		p_input_data++;
	}
	return rtn;
}

/**
 * @brief	Tests the next bits for zero groups
 * @param	p_num_bits		Pointer is set to the length of the next group of bits
 * @param	rle_limit		Minimum number of zeros before the stop condition is set
 * @return	1: if the next group is zero, 0: if the next group is not zero
 */
static uint8_t test_next_bits(uint_fast16_t* p_num_bits, uint8_t rle_limit)
{
	*p_num_bits = 0;

	uint8_t ref;
	uint8_t test_bits;
	uint8_t* _p_input_data;
	uint_fast8_t _input_bit_pos;
	uint8_t zero_cntr;

	ref = get_next_input_bits() == 0;
	do {
		_p_input_data = p_input_data;
		_input_bit_pos = input_bit_pos;
		test_bits = get_next_input_bits() == 0;
		(*p_num_bits)++;

		/* look ahead for zeros */
		zero_cntr = rle_limit - 1;
		while ((!ref && test_bits) && zero_cntr--){
			test_bits = get_next_input_bits() == 0;
			(*p_num_bits)++;
		}

	} while ( ref == test_bits && !end_of_file);

	zero_cntr++;
	if (!zero_cntr){
		*p_num_bits -= rle_limit - 1;
	}

	if (ref){
		p_input_data = _p_input_data;
		input_bit_pos = _input_bit_pos;
	}

	return ref;
}

/**
 * @brief	Writes a new zeros block to the rle data stream
 * @param	num_zeros		Number of zeros
 * @return	Nothing
 */
static void encode_zeros(uint16_t num_zeros)
{
	static uint16_t prev_num_zeros = 0;

	uint32_t new_data;
	uint_fast8_t new_data_len;

	uint16_t zero_blocks = num_zeros / 8;
	uint8_t zeros = num_zeros % 8;

	if ( (num_zeros - prev_num_zeros) > -7 && (num_zeros - prev_num_zeros) < 7){
		new_data = (num_zeros - prev_num_zeros) + 6;
		new_data_len = 1;
	} else if ( (num_zeros - prev_num_zeros) > -15 && (num_zeros - prev_num_zeros) < 15){
		new_data = 0xd << 4;
		if ((num_zeros - prev_num_zeros) >= 7 ){
			new_data |= (num_zeros - prev_num_zeros) + 1;
		} else {
			new_data |= (num_zeros - prev_num_zeros) + 14;
		}
		new_data_len = 2;
	} else if (zero_blocks < 32){
		new_data = 0xe << 8;
		new_data |= zero_blocks << 3 | zeros;
		new_data_len = 3;
	} else {
		new_data = 0xf << 16;
		new_data |= zero_blocks << 3 | zeros;
		new_data_len = 5;
	}

	uint8_t next_bits;
	while(new_data_len--){
		next_bits = (new_data >> (new_data_len * 4)) & 0xf;
		*p_rle_data |= next_bits << ((rle_bit_pos % 2) * 4);
		if (!(rle_bit_pos % 2)){
			p_rle_data++;
			*p_rle_data = 0;
		}
		rle_bit_pos++;
	}

	prev_num_zeros = num_zeros;
}

/**
 * @brief	Writes table encoded data to the stream
 * @param	mode		Compression mode
 * @param	num_bits	Size of the input data
 * @return	Nothing
 */
static void encode_data(enum comp_mode mode, uint16_t num_bits)
{
	static uint8_t val = 0;
	uint8_t next_val;
	uint8_t diff;
	uint8_t shifts;

	while(num_bits--){
		next_val = get_next_input_bits();
		switch (mode){
		case SUB_MODE:
			diff = (val - next_val) & 0xf;
			break;
		case MOV_MODE:
			diff = next_val & 0xf;
			break;
		case XOR_MODE:
			diff = (val ^ next_val) & 0xf;
			break;
		}
		shifts = compr_lut_shifts[diff];
		val = next_val;

		/* write the lut data to the current bit position */
		if (delta_bit_pos == 0){
			*p_delta_data = compr_lut_offset[diff];
		} else {
			*p_delta_data |= compr_lut_offset[diff] >> delta_bit_pos;
			*(p_delta_data + 1) |= compr_lut_offset[diff] << ( 8 - delta_bit_pos);
		}

		/* advance the pointer if the shift leaves the boundaries */
		delta_bit_pos += shifts;
		if (delta_bit_pos >= 8){
			delta_bit_pos -= 8;
			p_delta_data++;
		}
	}

	if (delta_bit_pos == 0){
		*p_delta_data = compr_offset_stop;
	} else {
		*p_delta_data |= compr_offset_stop >> delta_bit_pos;
		*(p_delta_data + 1) |= compr_offset_stop << ( 8 - delta_bit_pos);
	}

	delta_bit_pos += compr_shifts_stop;
	if (delta_bit_pos >= 8){
		delta_bit_pos -= 8;
		p_delta_data++;
	}

}

/**
 * @brief	Compresses images for Grandia
 * @param	mode			Compression mode to use
 * @param	input			Input file name as string
 * @param	output			Output file name as string
 * @param	p_table			Pointer to compression table
 * @param	rle_limit		RLE limit parameter
 * @return	Status
 */
compression_status_t compress(
		enum comp_mode mode,
		const char* input,
		const char* output,
		uint8_t* p_table,
		uint8_t rle_limit)
{
	int size = 0;
	int rsize = 0;
	FILE *f = fopen(input, "rb");

	uint8_t table[9];
	memcpy(table, p_table, 9);

	if (f == NULL){
		return FILE_NOT_FOUND;
	}

	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);

	p_input_data = (uint8_t *) malloc(size);
	p_end_of_file = p_input_data + size -1;
	rsize = fread(p_input_data, sizeof(uint8_t), size, f);

	if (size != rsize)
	{
		free(p_input_data);
		return INSUFFICIENT_MEM;
	}
	fclose(f);


	/* allocate memory for compressed data */
	uint8_t* p_delta_data_start;
	uint8_t* p_rle_data_start;
	p_delta_data_start = (uint8_t *) calloc(32000, sizeof(uint8_t*));
	p_rle_data_start = (uint8_t *) malloc(16000);
	p_delta_data = p_delta_data_start;
	p_rle_data = p_rle_data_start;

	/* generate look up table */
	if (generate_lut(p_table)){
		free(p_delta_data_start);
		free(p_rle_data_start);
		free(p_input_data);
		return INVALID_TABLE;
	}

	uint_fast16_t bit_cntr = 0;
	uint8_t is_zeros;
	*p_rle_data = 0;
	uint_fast32_t run_cntr = 0;
	uint16_t head;

	/* test if file starts on table encoded data */
	uint8_t initial_zeros;
	store_input_pos(0);
	is_zeros = test_next_bits(&bit_cntr, rle_limit);
	if (is_zeros && (bit_cntr < rle_limit)){
		initial_zeros = bit_cntr;
		test_next_bits(&bit_cntr, rle_limit);
		head = 0x2000;
		store_input_pos(1);
		encode_data(mode, bit_cntr+initial_zeros);
	} else {
		head = 0;
		store_input_pos(1);
	}

	while(!end_of_file){

		run_cntr++;
		if (run_cntr > MAX_COMP_LEN){
			free(p_delta_data_start);
			free(p_rle_data_start);
			free(p_input_data);
			return RUNAWAY_DATA;
		}

		store_input_pos(0);
		is_zeros = test_next_bits(&bit_cntr, rle_limit);

		if (is_zeros){
			encode_zeros(bit_cntr);
		} else {
			store_input_pos(1);
			encode_data(mode, bit_cntr);
		}
	}

	/* write remaining rle data and stop condition */
	uint32_t new_data = 0xe00;
	uint8_t new_data_len = 3;
	uint8_t next_bits;
	while(new_data_len--){
		next_bits = (new_data >> (new_data_len * 4)) & 0xf;
		*p_rle_data |= next_bits << ((rle_bit_pos % 2) * 4);
		if (!(rle_bit_pos % 2)){
			p_rle_data++;
			*p_rle_data = 0;
		}
		rle_bit_pos++;
	}
	if (!(rle_bit_pos % 2)){
		p_rle_data++;
	}

	/* write lut data stop condition if there is a trailing rle block */
	if (is_zeros){
		if (delta_bit_pos == 0){
			*p_delta_data = compr_offset_stop;
		} else {
			*p_delta_data |= compr_offset_stop >> delta_bit_pos;
			*(p_delta_data + 1) |= compr_offset_stop << ( 8 - delta_bit_pos);
		}

		delta_bit_pos += compr_shifts_stop;
		if (delta_bit_pos >= 8){
			delta_bit_pos -= 8;
			p_delta_data++;
		}
	}

	/* write remaining lut data */
	if (delta_bit_pos != 0){
		p_delta_data++;
	}

	/* write the header */
	switch (mode){
	case SUB_MODE:
		head |= 0b10 << 14;
		break;
	case MOV_MODE:
		head |= 0b00 << 14;
		break;
	case XOR_MODE:
		head |= 0b01 << 14;
		break;
	}
	head |= (uint16_t) (p_rle_data - p_rle_data_start);
	head = __bswap_16(head);

	f = fopen(output, "wb");
	fwrite(&head, sizeof(int16_t), 1, f);
	fwrite(p_rle_data_start, sizeof(int8_t), (int32_t) (p_rle_data - p_rle_data_start), f);
	fwrite(table, sizeof(int8_t), 9, f);
	fwrite(p_delta_data_start, sizeof(int8_t), (int32_t) (p_delta_data - p_delta_data_start), f);
	fclose(f);

	free(p_delta_data_start);
	free(p_rle_data_start);
	free(p_input_data);

	return SUCCESS;
}
