# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 16:48:05 2021

@author: nanashi
"""

import numpy as np
from . import lossless
from . import lossy
from . import colorspace

def encode_frame(frame : np.ndarray, level : int, subsample="discard", **kwargs):
    """
    Encodes a frame from an RGB numpy array format to the Grandia video stream
    format.

    Parameters
    ----------
    frame : np.ndarray
        Numpy array of the image frame. Shape must be (height, width, 3) with
        8 bit RGB values as entries.
    level : int
        Compression level for this frame. Must be an integer value 
        between 0 and 16. 0 is the lowest compression and 16 the highest.
    subsample : "discard" or "avrg"
        Determines how the chroma subsampling is handled. In discard mode
        only the upper left pixel in 4x4 unit is sampled. In average mode
        the upper left and lower left pixel are averaged together.
        
    Keywords
    ----------
    rclip : (int, int)
        Clips red channel colors to the specified boundaries before the colorspace
        conversion to avoid overflows.
    gclip : (int, int)
        Clips green channel colors to the specified boundaries before the colorspace
        conversion to avoid overflows.
    bclip : (int, int)
        Clips blue channel colors to the specified boundaries before the colorspace
        conversion to avoid overflows.
    yscale : (int, int)
        Scale y channel to the specified boundaries to avoid 
        overflows. y ranges from 0-255. See https://en.wikipedia.org/wiki/YCbCr
    cscale : int
        Scale cb and cr channel to the specified boundaries to avoid 
        overflows. C channels range from -128 to 127. Since c channels are 
        symmetric there is only one value.
        

    Returns
    -------
    bytearray
        Binary encoded data.

    """
    r = frame[:,:,0].astype(float)
    g = frame[:,:,1].astype(float)
    b = frame[:,:,2].astype(float)
    
    y, cb, cr = colorspace.rgb2ycbcr(r, g, b, subsample=subsample, **kwargs)
    
    y = lossless.gen_blocks(y)
    cr = lossless.gen_blocks(cr)
    cb = lossless.gen_blocks(cb)
    
    y = [[lossy.encode(block, level) for block in line] for line in y]
    cb = [[lossy.encode(block, level) for block in line] for line in cb]
    cr = [[lossy.encode(block, level) for block in line] for line in cr]
    
    return lossless.encode(y, cb, cr, level)