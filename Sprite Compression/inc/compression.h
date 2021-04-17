/*
 * compression.h
 *
 *  Created on: 15.09.2019
 *      Author: nanash1
 */

#ifndef INC_COMPRESSION_H_
#define INC_COMPRESSION_H_

#define MAX_DECOMP_LEN		(5000)
#define MAX_COMP_LEN		(10000)

extern const uint8_t lut_seed[9];

enum comp_mode {
	MOV_MODE,
	XOR_MODE,
	SUB_MODE,
};

typedef enum {
	SUCCESS,
	FILE_NOT_FOUND,
	INSUFFICIENT_MEM,
	RUNAWAY_DATA,
	INVALID_HEADER,
	INVALID_TABLE
} compression_status_t;

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
		uint8_t rle_limit);

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
		uint8_t* p_rle_limit);

#endif /* INC_COMPRESSION_H_ */
