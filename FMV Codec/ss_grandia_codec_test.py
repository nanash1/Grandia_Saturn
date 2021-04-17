# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 10:42:27 2021

@author: nanashi
"""
import numpy as np
from decord import VideoReader
from decord import cpu
from PIL import Image
import copy as cpy

from ss_grandia_codec import lossless
from ss_grandia_codec import lossy
from ss_grandia_codec import colorspace

'''
Settings
'''
sector = (6,11)
frame_no = 900
s1 = sector[0]*2
s2 = sector[1]*2
s3 = sector[0]*2+1
s4 = sector[1]*2+1

'''
Read decoded data from game dump of HWRAM
'''
decoded_y = np.zeros((176,352), dtype=int)
with open('test/esp_y.bin', 'rb') as decoded:
    for i in range(0, 176):
        for j in range(0, 352):
            decoded_y[i,j] = int.from_bytes(decoded.read(2), byteorder='big', signed=True)

decoded_cr = np.zeros((88,176), dtype=int)
with open('test/esp_cb.bin', 'rb') as decoded:
    for i in range(0, 88):
        for j in range(0, 176):
            decoded_cr[i,j] = int.from_bytes(decoded.read(4), byteorder='big', signed=True)

decoded_cb = np.zeros((88,176), dtype=int)            
with open('test/esp_cr.bin', 'rb') as decoded:
    for i in range(0, 88):
        for j in range(0, 176):
            decoded_cb[i,j] = int.from_bytes(decoded.read(4), byteorder='big', signed=True)

'''
Extract reference 16x16 sector
'''
game_cb = decoded_cb[sector[0]*8:(sector[0]*8)+8,sector[1]*8:(sector[1]*8)+8]
game_cr = decoded_cr[sector[0]*8:(sector[0]*8)+8,sector[1]*8:(sector[1]*8)+8]
game_y = np.block([[decoded_y[s1*8:(s1*8)+8,s2*8:(s2*8)+8], decoded_y[s1*8:(s1*8)+8,s4*8:(s4*8)+8]], 
                     [decoded_y[s3*8:(s3*8)+8,14*8:(14*8)+8], decoded_y[s3*8:(s3*8)+8,s4*8:(s4*8)+8]]])

'''
Read decoded data from game dump of VDP2RAM
'''
image_raw = bytes()
with open('test/esp_vdp2.bin', 'rb') as im:
    while True:
        data = im.read(4)
        
        image_raw += data[3:4]+data[2:3]+data[1:2]
        if not data:
            break

'''
Extract reference 16x16 sector
'''
game_rgb = np.empty((16,16), dtype='object')
for i in range(0,16):
    for j in range(0,16):
        game_rgb[i,j] = (image_raw[sector[0]*16896+(i*1056)+sector[1]*48+(j*3)], 
                           image_raw[sector[0]*16896+(i*1056)+sector[1]*48+(j*3)+1], 
                           image_raw[sector[0]*16896+(i*1056)+sector[1]*48+(j*3)+2])

'''
Read test frame #900 from video
'''
frame = VideoReader("examples/intro.mkv", ctx=cpu(0))[frame_no].asnumpy()
frame = Image.fromarray(frame, mode='RGB')
frame = frame.resize((352, 198), resample=Image.LANCZOS)
frame = frame.crop((0, 11, 352, 187))
frame = np.array(frame)
#imshow(frame)

'''
Covert RGB to YCbCr
'''
r = frame[:,:,0].astype(float)
g = frame[:,:,1].astype(float)
b = frame[:,:,2].astype(float)

y, cb, cr = colorspace.rgb2ycbcr(r, g, b, subsample="discard", rclip=(3,251), gclip=(3,251),bclip=(3,251))

'''
Generate macroblocks
'''
y = lossless.gen_blocks(y)
cr = lossless.gen_blocks(cr)
cb = lossless.gen_blocks(cb)

'''
Save source channels as reference
'''
src_y = np.block([[y[s1][s2], y[s1][s4]], 
                  [y[s3][s2], y[s3][s4]]])
src_cb = cpy.copy(cb[sector[0]][sector[1]])
src_cr = cpy.copy(cr[sector[0]][sector[1]])


'''
Apply lossy compression
'''
y = [[lossy.encode(block, 10) for block in line] for line in y]
cb = [[lossy.encode(block, 10) for block in line] for line in cb]
cr = [[lossy.encode(block, 10) for block in line] for line in cr]

'''
Decode the lossy compressed test sector
'''
dct_y0 = lossless.reduce_block(y[s1][s2])
dct_y1 = lossless.reduce_block(y[s1][s4])
dct_y2 = lossless.reduce_block(y[s3][s2])
dct_y3 = lossless.reduce_block(y[s3][s4])
dct_cb = lossless.reduce_block(cb[sector[0]][sector[1]])
dct_cr = lossless.reduce_block(cr[sector[0]][sector[1]])

dec_y = np.block([[lossy.decode(dct_y0, 10), lossy.decode(dct_y1, 10)], 
                  [lossy.decode(dct_y2, 10), lossy.decode(dct_y3, 10)]])
dec_cb = lossy.decode(dct_cb, 10)
dec_cr = lossy.decode(dct_cr, 10)

'''
Covert test sector to RGB
'''
test_rgb = np.empty((16,16), dtype='object')
for i in range(0,8):
    for j in range(0,8):
        
        if i == 7 and j == 6:
            test = 1
        res = colorspace.ycbcr2rgb_test(src_y[i*2,j*2], src_y[i*2,j*2+1], src_y[i*2+1,j*2], src_y[i*2+1,j*2+1], src_cb[i,j], src_cr[i,j])
        test_rgb[2*i,2*j] = res[0]
        test_rgb[2*i,2*j+1] = res[1]
        test_rgb[2*i+1,2*j] = res[2]
        test_rgb[2*i+1,2*j+1] = res[3]

'''
Covert test sector to RGB
'''
dec_rgb = np.empty((16,16), dtype='object')
for i in range(0,8):
    for j in range(0,8):
        
        res = colorspace.ycbcr2rgb_420(dec_y[i*2,j*2], dec_y[i*2,j*2+1], dec_y[i*2+1,j*2], dec_y[i*2+1,j*2+1], dec_cb[i,j], dec_cr[i,j])
        dec_rgb[2*i,2*j] = res[0]
        dec_rgb[2*i,2*j+1] = res[1]
        dec_rgb[2*i+1,2*j] = res[2]
        dec_rgb[2*i+1,2*j+1] = res[3]
        
'''
Test direct YCbCr to RGB conversion on game dump data
'''
game_dec_rgb = np.empty((16,16), dtype='object')
for i in range(0,8):
    for j in range(0,8):
        
        res = colorspace.ycbcr2rgb_420(game_y[i*2,j*2], game_y[i*2,j*2+1], game_y[i*2+1,j*2], game_y[i*2+1,j*2+1], game_cb[i,j], game_cr[i,j])
        game_dec_rgb[2*i,2*j] = res[0]
        game_dec_rgb[2*i,2*j+1] = res[1]
        game_dec_rgb[2*i+1,2*j] = res[2]
        game_dec_rgb[2*i+1,2*j+1] = res[3]

'''
Apply lossless copmression and write the generated frame data
'''
compressed = lossless.encode(y, cb, cr, 10)

# write frame
written = 0
with open('test.bin', 'wb') as outest:
    written += outest.write(compressed)
    if written % 0x800:
        padding = 0x800 - (written % 0x800)
        written += outest.write(bytearray(padding))