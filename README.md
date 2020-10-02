0. Disclaimer 

This software is supplied "AS IS" without any warranties of any kind. The author, nanash1, assumes no responsibility or liability for the use of the software, conveys no license or rights under any patent, copyright, mask work right, or any other intellectual property rights in or to any products. The author also makes no representation or warranty that such application will be suitable for the specified use without further testing or modification. 

Permission to use, copy, modify, and distribute this software and its documentation is hereby granted. This disclaimer notice must appear in all copies of this code. 

1. Introduction 

This tool is written to compress and decompress graphics for the original Saturn version of Grandia. 

The decompressed output images are encoded in a continuous binary file that contains the images as they are written to the Saturn's VDP1 Video RAM. This means if no color pallet is encoded in the compressed stream no pallet will be written to the output and all images within one stream are written to the same output file. In general, the images will not be viewable by conventional graphics editing software without conversion. A possible way to convert the output images is Magicks convert command. Since different kind images may be encoded in the same compressed stream, it's up to the user to separate them and find out properties like height, width, color depth, etc. before conversion. This is not within the scope of this program. 

Input images that are to be compressed must be laid out exactly like they were extracted. If they're not laid out in the same way, they will be written to VDP1's RAM incorrectly and appear as garbage in the game. It's also important to inject the newly compressed data at the right offset in the target file, so that Grandia begins to read the data at beginning. It should be noted that editing decompressed images and recompressing them may change the file size. There seems to be some padding in the CD data that may be used to expand image data, but at some point images may run beyond their boundaries and corrupt the game data. 

2. Usage

2.1 Decompression

To decompress images use the decompress mode of the program by using "decompress" as first argument: 

ss_grandia_comp_tool.exe decompress {-c [data offset]} -i [input file] -o [output file]
		
[input file]					Input file to compress/decompress.
[output file]					Output will be written to this file.
[compressed data offset]		Specifies where compressed data 
								begins within the input file. 
								Default: 0x0.

If the decompression was successful, the properties that were used for the decompression will be shown. This information is important for re-compressing the data. In general, Grandia should understand all valid compression options, however using inproper options will result in bad compression or even file enlargement. Unless the images were heavily edited it's probably best to trust the original developers and use the same options for the compression as for the decompression. 

Example for decompression:
ss_grandia_comp_tool.exe decompress -c 0xcd6b90 -i "Track 01.iso" -o image.bin

Result:
Output was written successfully.
Table: 0x145555566666545460
RLE limit: 6
MODE: SUB

2.2 Compression

To compress images use the compress mode of the program by using "compress" as first argument: 

ss_grandia_comp_tool.exe compress {[mode] -z [rle limit]} -t [table] -i [input file] -o [output file]

[input file]					Input file to compress/decompress.
[output file]					Output will be written to this file.
[table]							Decompression look up table.
[mode]							Defines the compression mode: 
									-s for SUB 
									-x for XOR 
									-m for MOV 
								Default: SUB.
[rle limit]						Minimum number of zeros in a RLE 
								encoded block. Default: 6.						

It's best to chose the same options for compression that were shown during the decompression. These options will usually yield the best compression possible. For more informations about the options please see the next section. 

Example for compression:
ss_grandia_comp_tool.exe compress -s -z 6 -t 0x145555566666545460 -i image.bin -o compressed.bin

Result:
Output was written successfully.
Table: 0x145555566666545460
RLE limit: 6
MODE: SUB

2.3 RLE limit

The RLE limit tells the compressor at which number of consecutive zeros a new RLE encoded block is started. Setting this option too low will result in too many RLE blocks and will increase the related overhead of the 4 bit wide identifiers. Setting it too high will result in less RLE encoded blocks and will encode more zeros in the look up table encoded stream. As this stream has in general less compression this also will result in a larger file size. 

2.4 Compression mode

Data values (each 4 bits wide) of the image can be written in three different ways to the look up table encoded stream: They can be written directly to the table. This is called MOV mode. The difference to the previous value can be written. This is called SUB mode. Or the value can be XOR'd with the previous value. XOR mode. What mode yields the best compression depends on the input image. For example, smooth gradients will probably compress best with SUB mode. 

2.5 Table

The compression table maps every possible value of a 4 bit integer to a bit length. Looking at the example table "0x145555566666545460" the way it works is this: The first value is 0x1, this means 0 is mapped to 1 bit. The second value is 0x4, this means 1 is mapped to 4 bits. The third value is 0x5, this means 2 is mapped to 5 bits. And so on. When Grandia decompresses the image it will generate an array with a size of 256 byte sized elements. Each number of bytes requires a specific number of elements to make each bit value identifiable in a 8 bit sliding window. 1 bit requires 128 elements, 2 bits require 64 elements, 3 bits require 32 elements and so on. Each table is 16 times 4 bits long plus 4 bits for the stop condition and plus 4 bits of zeros for memory alignment. If a value doesn't appear in the the data, it can be be skipped in the table by setting the bit length to 0.

This means the following requirements for a valid table:
- 9 bytes long, 4 LSBs are zero
- Results in 256 elements after expansion
- Maximum number of bits is 8

Invalid example
0x14555556666654547
	- table not 256 elements long
	- zero padding is missing

Valid examples
0x145555566666545460
0x555555555344444250
0x043334500445444440

To get a good compression it's important to assign frequent valuesto a low number of bits. For example, the table from the exampleassigns 0x0 to 1 bit. That's because the table is optimized for SUB mode compression and a difference of zero is most common. 

3. Known offsets within extracted user data of the game image

Disc 1:
0x25ebc				Battle popup images
0xfb3a48			Battle results screen images
0x1003218			Battle wheel images

Disc 2:
0x266bc				Battle popup images
0x1a26a48			Battle results screen images
0x1a76218			Battle wheel images

Digital Museum Disc:
0x24634				Battle popup images
0xc87248			Battle results screen images
0xcd6b90			Battle wheel images

4. Brief description of the compression used

Grandia uses a combination of RLE and a look up table to compress images. First RLE is used to compress all zeros and then the look up table is used on the remaining data. The outputs of these compressions are stored separately in the compressed data. The layout is as follows

1. 16 bit header:			Contains delta encoding type (XOR, SUB or MOV) and the offset at which the
							table encoded data begins relative to the header
2. RLE data:				Compressed zeros data in blocks
3. Decompression table:		Contains mappings of values to bits (see above)
4. Table encoded data:		Table encoded data stripped of zeros
	
4.1 RLE data

Example data: E81 E24 8 8 8 8

As stated before the RLE data contains the encoded zeros. The encoded zero blocks can be either 16 bits, 8 bits and 4 bits. They can either be absolute or differential. The type of the block is specified by a prefix contained in the first 4 bits. 0xF means 16 bits absolute, 0xE means 8 bits absolute, 0xD means 4 bits differential with offset and no prefix means 4 bits differential without offset. Differential means only the difference in zeros compared to the last block is stored. Let's look at the example data above and decode it:

The first block is E81. This means there are 0x81 zeros or 129. The next block is E24, so there are 36 zeros. The next blocks don't have a prefix so they are differential. The value is 0x8, but to get the offset, 6 is subtracted to allow a negative difference. So this block has 36 + (8 - 6) = 38 zeros. The next value is also 0x8. This means the block has 38 + (8 - 6) = 40 zeros and so forth.

4.2 Table encoding (example in SUB mode)

Example image data:		6EF888888888888888875	
Example encoded data: 	0C1703FFFA8C

The table encoded data uses a look up table that is compressed in memory to map deltas to bits. Each value is read as 8 bits from the encoded stream. This means values are overlapping in the encoded stream depending on their bit size. To keep values identifiable, the expanded table needs an increasing number of elements. For example if a value is encoded in 1 bit, it has 7 remaining bits that can beany value. This means 2^7 = 128 elements are needed in the expanded table.

1. 0-6 = A; from LUT 0xA := 0b0000|11
	=> Output 0b0000|11
2. 6-E = 8; from LUT 0x8 := 0b0000|01
	=> Output 0b0000|1100|0001|
3. E-F = F; from LUT 0xF := 0b0111
	=> Output 0b0000|1100|0001|0111|
4. F-8 = 7; from LUT 0x7 := 0b0000|00 
	=> Output 0b0000|1100|0001|0111|0000|00
5. 8-8 = 0; from LUT 0x0 := 0b1
	=> Output 0b0000|1100|0001|0111|0000|001
6. 8-8 = 0; from LUT 0x0 := 0b1
	=> Output 0b0000|1100|0001|0111|0000|0011|
7. 8-8 = 0; from LUT 0x0 := 0b1
	=> Output 0b0000|1100|0001|0111|0000|0011|1
	
And so forth.

5. Version history

v0.95
	- Fixed unhandled case when images start on table encoded data
v0.9
	- Source code cleanup
v0.8
	- Fixed XOR decompression
v0.7
	- Added handling of look up tables with unmapped values
	- Improved look up table validity check
	- Fixed bug if RLE data starts on an offset block

v0.6
	- Added support for MOV and XOR mode
	- Hardcoded table seed data since it's always the same
	- Added custom compression tables to compression function
	- Added selectable mode for compression
	- Added selectable RLE limit for compression

v0.5 Initial release
	- XOR compression not implemented yet
	- Custom compression look up tables not supported for compression yet
