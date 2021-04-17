/*
 ============================================================================
 Name        : ss_grandia_comp_tool.c
 Author      : nanash1
 Version     : 0.95
 Description : Tool to compress and decompress Grandia Saturn graphics
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
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "compression.h"

enum mode {
		COMPRESS,
		DECOMPRESS
};

uint8_t lut_table[9];
char lut_table_str[22] = {0};

static void lut_to_string(void)
{
	lut_table_str[0] = '0';
	lut_table_str[1] = 'x';
	for (uint_fast8_t i = 0; i < 9; i++){
		if (lut_table[i] < 0xf){
			lut_table_str[2+2*i] = '0';
			sprintf(&lut_table_str[2] + 2*i+1, "%x", lut_table[i]);
		} else {
			sprintf(&lut_table_str[2] + 2*i, "%x", lut_table[i]);
		}
	}
}

int main(int argc, char* argv[])
{
	compression_status_t status;
	enum mode op_mode;
	uint_fast8_t arg_found;
	const char* input;
	const char* output;

	uint32_t compressed_data_offset = 0x0;
	enum comp_mode comp_mode = SUB_MODE;
	uint8_t rle_limit = 6;

	uint8_t table_arg_ind;
	size_t table_size;

	printf("Saturn Grandia graphics compression tool version v0.95 by nanash1\n\n");

    if(argc == 1){
    	printf("No arguments specified. Please refer to the readme file.\n");
    	return 0;
    }

    if(!strcmp(argv[1], "compress")){
    	op_mode = COMPRESS;
    } else if (!strcmp(argv[1], "decompress")){
    	op_mode = DECOMPRESS;
    } else {
    	printf("Invalid operations mode argument. Mode must be compress or decompress.\n");
    	return 0;
    }

    arg_found = 0;
    for (uint_fast8_t i = 1; i < (argc - 1); i++){
    	if (!strcmp(argv[i], "-i")){
    		input = argv[i+1];
    		arg_found++;
    		break;
    	}
    }
    for (uint_fast8_t i = 1; i < (argc - 1); i++){
    	if (!strcmp(argv[i], "-o")){
    		output = argv[i+1];
    		arg_found++;
    		break;
    	}
    }

    if (arg_found != 2){
    	printf("Please specify input and output file.\n");
    	return 0;
    }

    if (op_mode == DECOMPRESS){

        arg_found = 0;
        for (uint_fast8_t i = 1; i < (argc - 1); i++){
        	if (!strcmp(argv[i], "-c")){
        		compressed_data_offset = strtol(argv[i+1], NULL, 0);
        		break;
        	}
        }

    	if (errno){
    		printf("Error: %d.\n", errno);
    		return 0;
    	}

        status = decompress(input, output, compressed_data_offset, &comp_mode, &lut_table[0], &rle_limit);

    } else if (op_mode == COMPRESS){

        arg_found = 0;
        for (uint_fast8_t i = 1; i < (argc - 1); i++){
        	if (!strcmp(argv[i], "-t")){

        		table_arg_ind = i + 1;
        		/* check if the table has a valid format */
        		table_size = strlen(argv[i+1]);
        		char segment[3] = {0};
        		memcpy(&segment[0], argv[i+1], 2);
        		if (strcmp(segment, "0x")){
                	printf("Table must be specified as hex value. Use the prefix 0x.\n");
                	return 0;
        		}
        		if (table_size < 20){
                	printf("Table must be at least 9 bytes long. Note that further bytes are ignored.\n");
                	return 0;
        		}

        		/* copy the table data to memory */
        		for (uint_fast8_t j = 0; j < 9; j++){
        			memcpy(&segment[0], &argv[i+1][2+(2*j)], 2);
        			lut_table[j] = strtol(&segment[0], NULL, 16);
        		}
        		arg_found++;
        		break;
        	}
        }

        if (!arg_found){
        	printf("Please specify a decompression table.\n");
        	return 0;
        }

        for (uint_fast8_t i = 1; i < (argc - 1); i++){
        	if (!strcmp(argv[i], "-x")){
        		comp_mode = XOR_MODE;
        		break;
        	} else if (!strcmp(argv[i], "-s")){
        		comp_mode = SUB_MODE;
        		break;
        	} else if (!strcmp(argv[i], "-m")){
        		comp_mode = MOV_MODE;
        		break;
        	}
        }

        for (uint_fast8_t i = 1; i < (argc - 1); i++){
        	if (!strcmp(argv[i], "-z")){
        		rle_limit = (uint8_t) strtol(argv[i+1], NULL, 0);;
        		break;
        	}
        }

    	if (errno){
    		printf("Error: %d.\n", errno);
    		return 0;
    	}

        status = compress(comp_mode, input, output, (uint8_t*) &lut_table[0], rle_limit);
    }

    /* print informations if successful */
    if (status == SUCCESS){
    	printf("Output was written successfully.\n");

		if (op_mode == DECOMPRESS){
			lut_to_string();
			printf("Table: %s\n", lut_table_str);
		} else {
			printf("Table: %s\n", argv[table_arg_ind]);
		}
		printf("RLE limit: %u\n", rle_limit);
		printf("Mode: ");
		switch (comp_mode){
		case XOR_MODE:
			printf("XOR\n");
			break;
		case MOV_MODE:
			printf("MOV\n");
			break;
		case SUB_MODE:
			printf("SUB\n");
			break;
		}
	}

    /* print error messages */
    switch (status){
    case SUCCESS:
    	break;
    case FILE_NOT_FOUND:
    	printf("File not found.\n");
    	break;
    case INSUFFICIENT_MEM:
    	printf("Out of memory.\n");
    	break;
    case RUNAWAY_DATA:
    	printf("Couldn't find end of stream. Please check offsets or image data.\n");
    	break;
    case INVALID_HEADER:
    	printf("Invalid header. Please check the offsets.\n");
    	break;
    case INVALID_TABLE:
    	printf("Compression table invalid.\n");
    	break;
    }

    return 0;
}
