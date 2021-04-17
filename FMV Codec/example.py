# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 18:14:14 2021

@author: nanashi
"""

import numpy as np
import ss_grandia_codec
from decord import VideoReader
from decord import cpu
from PIL import Image

compression_level = 10
frame_number = 900

# Import video with decord's VideoReader
frame = VideoReader("examples/intro.mkv", ctx=cpu(0))[frame_number].asnumpy()

# Read into PIL.Image for downscaling and croping
frame = Image.fromarray(frame, mode='RGB')
frame = frame.resize((352, 198), resample=Image.LANCZOS)
frame = frame.crop((0, 11, 352, 187))

# Convert to numpy RGB array
frame = np.array(frame)

# Encode frame to raw video data
frame_data = ss_grandia_codec.encode_frame(frame, compression_level, subsample="avrg", yscale=(0,220), cscale=220)

# Write the result to a file
with open("encoded.bin", "wb") as outf:
    outf.write(frame_data)