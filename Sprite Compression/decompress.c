/*
 ============================================================================
 Name        : decompress.c
 Author      : nanash1
 Version     : 0.95
 Description : Decompresses Grandia Saturn graphics files
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

const uint8_t lut_seed[9] = {0, 128, 64, 32, 16, 8, 4, 2, 1};
static uint16_t decomp_lut[256];

/**
 * @brief  	Gets the next byte of encoded data and advances the pointer
 * @param	p_offset_code		Next bit offset to advance the pointer
 * @param	pp_lut_data			Pointer to current input data pointer
 * @return	Encoded data
 */
static uint16_t extract_lut_data(uint16_t* p_offset_code, uint8_t** pp_lut_data)
{
	uint16_t shift_code = *p_offset_code >> 4;
	uint32_t fetched_source_data;

	if (shift_code > 31){
		*pp_lut_data += 4;
		shift_code -= 32;
		*p_offset_code = shift_code << 4;
	}

	uint8_t* p_source_data = *pp_lut_data;

	if (shift_code < 25){
		fetched_source_data = __bswap_32(*((uint32_t*) p_source_data)) >> (24 - shift_code);
	} else {
		fetched_source_data = __bswap_32(*((uint32_t*) p_source_data)) << (shift_code - 24);
		fetched_source_data |= __bswap_32(*((uint32_t*) (p_source_data + 4))) >> (56 - shift_code);	// or in bits from next pointer
	}

	return (uint16_t) ((fetched_source_data & 0xff) << 1);
}

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
	uint16_t lut_elem;
	uint16_t *p_decomp_lut = &decomp_lut[0];

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
			break;
		}

		if (decode_char%2){
			*(p_comp_lut + (decode_char/2)) &= 0xf0;
		} else {
			*(p_comp_lut + (decode_char/2)) &= 0xf;
		}

		lut_elem = decode_char << 8;
		lut_elem |= max_bit_num << 4;
		if (decode_char == 16){
			lut_elem++;
		}
		for (uint_fast8_t i = 0; i < lut_seed[max_bit_num]; i++){

			*p_decomp_lut = lut_elem;
			p_decomp_lut++;
		}
	}

	/* check if the table is within bounds */
	if (p_decomp_lut != &decomp_lut[0] + 256){
		rtn = 1;
	}
	return rtn;

}

/**
 * @brief	Decompresses the compressed data from the Grandia iso
 * @param	input					Input file name as string
 * @param	output					Output file name as string
 * @param	comp_data_offset		Start offset of the compressed data in the input file
 * @param	p_mode					The used decompression mode is written to this pointer
 * @param	p_lut					The used look up table is written to this pointer
 * @param	p_rle_limit				The used RLE limit is written to this pointer
 * @return	Status
 */
compression_status_t decompress(
		const char* input,
		const char* output,
		uint32_t comp_data_offset,
		enum comp_mode* p_mode,
		uint8_t* p_lut,
		uint8_t* p_rle_limit)
{
	uint8_t* p_iso;


	/* read dumped file into memory */
	int size = 0;
	FILE *f = fopen(input, "rb");

	if (f == NULL){
		return FILE_NOT_FOUND;
	}

	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);

	p_iso = (uint8_t *) malloc(size+1);

	if (size != fread(p_iso, sizeof(int8_t), size, f))
	{
		free(p_iso);
		return INSUFFICIENT_MEM;
	}
	fclose(f);

	/* read header */
	uint8_t* p_rle_data = p_iso + comp_data_offset;
	uint16_t head_data;
	head_data = __bswap_16(*((uint16_t*) p_rle_data));

	/* get mode */
	uint8_t mode_bits = head_data >> 14;
	enum comp_mode mode;
	switch (mode_bits){
	case 0b10:
		mode = SUB_MODE;
		break;
	case 0b00:
		mode = MOV_MODE;
		break;
	case 0b01:
		mode = XOR_MODE;
		break;
	default:
		free(p_iso);
		return INVALID_HEADER;
	}
	*p_mode = mode;

	/* write decompression table */
	uint16_t table_offset = head_data;
	table_offset = table_offset & 0x1fff;
	uint8_t* p_table_data = p_rle_data + table_offset + 2;
	memcpy(p_lut, p_table_data, 9);
	if (generate_lut(p_table_data)){
		free(p_iso);
		return INVALID_TABLE;
	}

	/* reserve memory for decompressed data */
	uint8_t* p_decompressed = (uint8_t*) calloc(64000, sizeof(uint8_t));
	uint8_t* p_output = p_decompressed;

	/* get pointer to lut encoded data */
	uint8_t* p_lut_data = p_rle_data + table_offset + 0xb;
	p_lut_data = (uint8_t*) (((uint32_t) p_lut_data) & 0xfffffffc);	// probably to align for 32bit memory access, TODO: check if this works here

	/* advance pointer to first RLE data */
	p_rle_data += 2;

	/* initialize decompression vars */
	uint8_t rle_zeros_cntr = 0;
	uint8_t rle_limit = 0;
	uint16_t zero_num = 0;
	uint16_t zero_diff;
	uint32_t decoded_bits_num = 0;

	uint16_t shift_code = ((table_offset + 0xb) & 0x3) << 7;
	uint16_t decompressed_bits = 0;
	uint8_t rle_prefix;
	uint32_t lut_index = 0;
	uint16_t lut_data;
	uint_fast8_t is_aligned = 0;

	uint_fast32_t run_cntr = 0;

	if ( (head_data >> 8) & 0x20){		// if the second bit is set
		goto skip_first_zeros;
	}

	while(1){

		run_cntr++;

		if (run_cntr > MAX_DECOMP_LEN){
			free(p_decompressed);
			free(p_iso);
			return RUNAWAY_DATA;
		}

		/* get next rle bits */
		if (is_aligned){
			rle_prefix = *p_rle_data & 0xf;
			p_rle_data++;
			is_aligned = 0;
		} else {
			rle_prefix = *p_rle_data >> 4;
			is_aligned = 1;
		}

		/* check the prefix code */
		switch (rle_prefix){
			case 0xD:	// read next 4 bits as extended delta information
				if (is_aligned){
					zero_diff = *p_rle_data & 0xf;
					is_aligned = 0;
					p_rle_data++;
				} else {
					zero_diff = *p_rle_data >> 4;
					is_aligned = 1;
				}
				if (zero_diff >= 8){
					zero_num += (zero_diff - 1);
				} else {
					zero_num += zero_diff - 14;
				}
				break;

			case 0xE:	// read next 8 bits as number of zeros
				if (is_aligned){
					zero_num = (*p_rle_data & 0xf) << 4;
					p_rle_data++;
					zero_num |= *p_rle_data  >> 4;
				} else {
					zero_num = *p_rle_data;
					p_rle_data++;
				}
				break;

			case 0xF:	// read next 16 bits as number of zeros
				if (is_aligned){
					zero_num = (*p_rle_data & 0xf) << 12;
					p_rle_data++;
					zero_num |= *p_rle_data << 4;
					p_rle_data++;
					zero_num |= (*p_rle_data & 0xf0) >> 4;
				} else {
					zero_num = *p_rle_data << 8;
					p_rle_data++;
					zero_num |= *p_rle_data;
					p_rle_data++;
				}
				break;

			default:	// read prefix as number of zeros offset
				zero_num += rle_prefix - 6;
				break;
		}


		/* write the determined number of zeros */
		if (zero_num){

			if (decoded_bits_num % 2){
				p_output += (zero_num - 1) / 2;
			} else {
				p_output += zero_num / 2;
			}
			decoded_bits_num += zero_num;

		} else {
			break;
		}

		while (1){

skip_first_zeros:
			/* use the generated table to read actual data from memory */
			lut_index = extract_lut_data(&shift_code, &p_lut_data);
			lut_data = decomp_lut[lut_index/2];

			/* look for stop condition and write data to output */
			if (!(lut_data & 0x1)){
				shift_code += lut_data & 0xff;
				lut_data >>= 8;
				switch (mode){
				case SUB_MODE:
					decompressed_bits = (decompressed_bits - lut_data) & 0xf;
					break;
				case MOV_MODE:
					decompressed_bits = lut_data & 0xf;
					break;
				case XOR_MODE:
					decompressed_bits = (decompressed_bits ^ lut_data) & 0xf;
					break;
				}

				if (decompressed_bits == 0){
					rle_zeros_cntr++;
				} else {
					if (rle_limit < rle_zeros_cntr){
						rle_limit = rle_zeros_cntr;
					}
					rle_zeros_cntr = 0;
				}

				if (decoded_bits_num % 2){
					*p_output |= decompressed_bits;
					p_output++;
					decoded_bits_num++;
				} else {
					*p_output = decompressed_bits << 4;
					decoded_bits_num++;
				}

			} else {
				if (decoded_bits_num % 2){
					p_output++;
				}
				shift_code += lut_data & 0xf0;
				break;
			}
		}
	}

	f = fopen(output, "wb");
	fwrite(p_decompressed, sizeof(uint8_t), (int32_t) (p_output - p_decompressed), f);
	fclose(f);

	free(p_decompressed);
	free(p_iso);

	*p_rle_limit = rle_limit + 1;

	return SUCCESS;
}
