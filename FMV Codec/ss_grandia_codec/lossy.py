# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 14:29:45 2021

@author: nanashi
"""

import numpy as np
import importlib.resources
from . import dct

quant_file = importlib.resources.open_binary("ss_grandia_codec","quantize.npy")
quant = np.load(quant_file)
    
def encode(block, level : int):
    '''
    Executes lossy block encoding:
        - DCT
        - Quantization

    Parameters
    ----------
    block : np.ndarray
        Macroblock matrix.
    level : int
        Compression level
        

    Returns
    -------
    np.ndarray
        Processed block.

    '''
    block = dct.dct_f(block)
    block /= quant[level]
    block = np.round(block).astype(int)
    return block

def decode(block, level : int):
    '''
    Decodes lossy encoded block

    Parameters
    ----------
    block : np.ndarray
        Macroblock matrix.
    level : int
        Compression level.

    Returns
    -------
    np.ndarray
        Decoded block.

    '''
    
    block = block*quant[level][:block.shape[0],:block.shape[1]].astype('int32')
    return dct.idct_i(block)