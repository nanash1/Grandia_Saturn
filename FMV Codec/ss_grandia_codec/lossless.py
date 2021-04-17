# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 14:52:24 2021

@author: nanashi
"""

import numpy as np
import math
import importlib.resources

order_file = importlib.resources.open_binary("ss_grandia_codec","order.npy")
order = np.load(order_file, allow_pickle=True).reshape((8,8))

def gen_blocks(chan, bsize=8):
    """
    Generates macro blocks from from channel data

    Parameters
    ----------
    chan : Array of float
        Channel data as numpy array.
    bsize : int, optional
        Block size. The default is 8.

    Returns
    -------
    res : list
        Nested list with blocks.

    """
    block_width = int(chan.shape[1]/bsize)
    block_height = int(chan.shape[0]/bsize)
    
    lines = np.split(chan, block_height)
    blocks = []
    for i in range(0, block_height):
        blocks.append(np.split(lines[i], block_width, axis=1))
    return blocks

def reduce_block(block):
    '''
    Reduces the block size by finding the non-zero subblock of a block

    Parameters
    ----------
    block : numpy array
        Macroblock matrix.

    Returns
    -------
    numpy array
        Macroblock matrix.

    '''
    x, y = np.nonzero(block)
    if len(x) > 0:
        x = x.max()
    else:
        x = 0
    if len(y) > 0:
        y = y.max()
    else:
        y = 0
    return block[:x+1,:y+1]  

def gen_code(num_zeros, elem, idx, codes, nums):
    '''
    Generates block compression codes

    Parameters
    ----------
    num_zeros : int
        Number of unwritten zeros.
    elem : int
        Element value.
    idx : int
        index of element in flattened macroblock.
    codes : list
        new codes are appended to this list.
    nums : list of (int, int)
        Additional data is appended to this list (number, size in bits).

    Returns
    -------
    None.

    '''
    abs_elem = abs(elem)
    
    full = 0
    part = 0
    if num_zeros > 0:
        full = int(num_zeros / 9)
        part = num_zeros % 9
        
    codes += [29]*full
    nums += [(7,3)]*full                                                        # offset by 2, so 7 := 9 zeros
    
    if abs_elem == 1:
        if part == 5:
            codes += [28]
            part = 4
        elif part > 5:
            codes += [29]
            nums += [(part-2,3)]
            part = 0
            
        if elem == -1:
            part += 5
        codes += [part]
        return
        
    elif abs_elem == 2:
        if part > 2:
            codes += [29]
            nums += [(part-2,3)]
            part = 0
        elif part == 2:
            codes += [28]
            part = 1
            
        code = 10
        if elem == -2:
            code += 1
        if part == 1:
            code += 2
        codes += [code]
        return
        
    elif part > 0:
        if part == 1:
            codes += [28]
        else:
            codes += [29]
            nums += [(part-2,3)]
            
    if abs_elem < 5:
        code = 14
        if abs_elem == 4:
            code += 2
        
    elif abs_elem < 7:
        code = 18
        nums += [(abs_elem-5, 1)]
        
    elif abs_elem < 11:
        code = 20
        nums += [(abs_elem-7, 2)]
        
    elif abs_elem < 19:
        code = 22
        nums += [(abs_elem-11, 3)]
        
    elif abs_elem < 35:
        code = 24
        nums += [(abs_elem-19, 4)]
    else:
        code = 26
        if idx == 0:
            nums += [(abs_elem-35, 11)]
        elif idx < 3:
            nums += [(abs_elem-35, 10)]
        else:
            nums += [(abs_elem-35, 8)]
            
    if elem < 0:
        code += 1 
    codes += [code]

def gen_comp_data(block, codes, nums, dim_codes, last_elem):
    '''
    Generates compressed data from blocks

    Parameters
    ----------
    block : numpy array
        preproccessed macro block.
    codes : list of int
        Compression codes are appended to this list.
    nums : list of (int, int)
        Additional data is appended to this list (number, size in bits).
    dim_codes : list of int
        Dimension codes are appended to this list.
    last_elem : int
        previous first element.

    Returns
    -------
    first_elem : int
        first element.

    '''
    block = reduce_block(block)
    rows = block.shape[0]
    cols = block.shape[1]
    length = cols*rows
    dim_codes += [(cols-1 << 3) | (rows-1)]
    block = block.T.reshape(length)
    block = block[order[cols-1, rows-1]]
    zero_cntr = 0
    first_elem = block[0]
    block[0] -= last_elem
    for i in range(0,length):
        elem = block[i]
        if elem == 0:
            zero_cntr += 1
        else:
            gen_code(zero_cntr, elem, i, codes, nums)
            zero_cntr = 0
    codes += [30]
    return first_elem

def _gen_start_cond(n:int):
    '''
    Generates the starting weights that are ideal if all elements are equally
    distributed

    Parameters
    ----------
    n : int
        Number of elements.

    Returns
    -------
    np.array
        Array with the weight of each element.

    '''
    x = int(math.log(n)/math.log(2))
    i = n - 2**x
    
    return np.array((n-2*i)*[1/(2**x)]+(2*i)*[1/(2**(x+1))])

def _find_weights(counts):
    '''
    Finds the ideal weights for each element so that the compression is minimal

    Parameters
    ----------
    counts : dict
        Dictonary that contains the elements as keys and their number of
        occurrence as value.

    Returns
    -------
    codes : numpy array
        Codes as numpy array.
    weights : numpy array
        Weights for each code as numpy array.

    '''
    codes = np.array(list(counts.keys()))
    counts = np.array(list(counts.values()))
    sort_idx = np.argsort(counts)[::-1]
    codes = codes[sort_idx]
    counts = counts[sort_idx]
    weights = _gen_start_cond(len(codes))
    
    log2 = math.log(2)
    for i in range(0, len(counts)-1):
        pool = []
        cref = counts[i]
        wsum = 0
        csum = 0
        min_weight = 1
        for j in range(len(counts)-1, i, -1):
            if weights[j] > 1/256:
                if csum + counts[j] > cref:
                    break
                else:
                    csum += counts[j]
                    wsum += weights[j]
                    pool.append(j)
                    if weights[j] < min_weight:
                        min_weight = weights[j]
        if wsum == 0:
            break
        mult = wsum / weights[i]
        mult = int(math.log(mult)/log2)
        if mult < 1:
            continue
        mult = 2**mult
        while (min_weight / mult) < (1/256):
            mult /= 2
        wref = weights[i]*mult
        wsum = 0
        cntr = 0
        for idx in pool:
            cntr += 1
            wsum += weights[idx]
            if wsum == wref:
                weights[i] *= mult
                weights[pool[:cntr]] /= mult
                break
    
    weights = np.log(1/weights)/log2
    return codes, weights

def gen_sec0(sec2_linepos, sec3_linepos):
    sec0 = bytearray(66)
    sec0[1] = 16
    
    bytepos = 6
    for i in range(0, len(sec2_linepos)-1):
        sec0[bytepos] = (sec2_linepos[i][0] + 16) >> 8                          # +16 because the 16 bytes of the decoding table are included
        sec0[bytepos+1] = (sec2_linepos[i][0] + 16) & 0xff
        sec0[bytepos+2] = (sec3_linepos[i][0]) >> 8
        sec0[bytepos+3] = (sec3_linepos[i][0]) & 0xff
        sec0[bytepos+4] = sec2_linepos[i][1]
        sec0[bytepos+5] = sec3_linepos[i][1]
        
        bytepos += 6
    return sec0

def gen_sec1(y0, y1, cb, cr):
    sec1 = bytearray(1452)
    bytepos = 0
    for i in range(0, len(cb)):
        sec1[bytepos] = y0[2*i] << 2
        sec1[bytepos] |= y0[2*i+1] >> 4
        sec1[bytepos+1] = (y0[2*i+1] << 4) & 0xff
        sec1[bytepos+1] |= y1[2*i] >> 2
        sec1[bytepos+2] = (y1[2*i] << 6) & 0xff
        sec1[bytepos+2] |= y1[2*i+1]
        sec1[bytepos+3] = cr[i] << 2
        sec1[bytepos+3] |= cb[i] >> 4
        sec1[bytepos+4] = (cb[i] << 4) & 0xff
        
        bytepos += 6 
    return sec1

def gen_sec2(code_lines):
    
    # analyse the generated codes to find best compression
    code_count = {}
    for line in code_lines:
        for code in line:
            if code not in code_count:
                code_count[code] = 1
            else:
                code_count[code] += 1
    codes, weights = _find_weights(code_count)
    
    sec2_codes = bytearray(16)
    sec2 = bytearray(1048576)
    
    # generate binary code lut
    sort_idx = np.argsort(codes)
    codes = codes[sort_idx]
    weights = weights[sort_idx]
    bin_lut = {}
    addr = 0
    lut_repeats = (0, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01)
    while True:
        idx = np.argmax(weights)
        max_code = int(weights[idx])
        if max_code == 0:
            break
        code = codes[idx]
        if code % 2:
            sec2_codes[int(code/2)] |= int(weights[idx])
        else:
            sec2_codes[int(code/2)] |= int(weights[idx]) << 4
        bin_lut[code] = (addr >> (8-max_code), max_code)
        addr += lut_repeats[max_code]
        weights[idx] = 0
        
    bitpos= 0
    bytepos = 0
    line_pos = []
    for line in code_lines:
        for code in line:
            shift = 8 - bitpos - bin_lut[code][1]
            bitpos += bin_lut[code][1]
            if shift < 0:
                sec2[bytepos] |= bin_lut[code][0] >> -shift
                bytepos += 1
                bitpos = -shift
                shift = 8 + shift
            sec2[bytepos] |= (bin_lut[code][0] << shift) & 0xff
            
        if bitpos == 8:
            line_pos.append((bytepos+1, 0))
        else:
            line_pos.append((bytepos, bitpos))
        
    return sec2_codes+sec2[:bytepos+1], line_pos
    

def gen_sec3(num_lines):
    sec3 = bytearray(1048576)
    bitpos= 8
    bytepos = 0
    line_pos = []
    for line in num_lines:
        for pair in line:
            bit_size = pair[1]
            while bit_size:
                shift = bit_size - bitpos
                if shift < 0:
                    sec3[bytepos] |= (pair[0] << -shift) & 0xff
                    bits_written = bit_size
                else:
                    sec3[bytepos] |= (pair[0] >> shift) & 0xff
                    bits_written = bitpos
                    bytepos += 1
                bitpos -= bits_written
                bit_size -= bits_written
                if bitpos == 0:
                    bitpos = 8
                    
        line_pos.append((bytepos, 8-bitpos))
    return sec3[:bytepos+1], line_pos

def encode(y, cb, cr, level):
    
    width = int(len(y[0])/2)
    height = int(len(y)/2)
    comp_codes = []
    comp_nums = []
    
    dim_codes_y0 = []
    dim_codes_y1 = []
    dim_codes_cb = []
    dim_codes_cr = []
    
    # generate compression codes and data for each line of the video frame
    # store dimension codes seperately because they aren't stored by line
    for i in range(0, len(cb)):
        comp_code_line = []
        comp_num_line = []
        diff = 0
        y_ind = int(i*2)
        for block in y[y_ind]:
            diff = gen_comp_data(block, comp_code_line, comp_num_line, dim_codes_y0, diff)
        diff = 0
        for block in y[y_ind+1]:
            diff = gen_comp_data(block, comp_code_line, comp_num_line, dim_codes_y1, diff)
        diff = 0
        for block in cr[i]:
            diff = gen_comp_data(block, comp_code_line, comp_num_line, dim_codes_cr, diff)
        diff = 0
        for block in cb[i]:
            diff = gen_comp_data(block, comp_code_line, comp_num_line, dim_codes_cb, diff)
        comp_codes.append(comp_code_line)
        comp_nums.append(comp_num_line)
        
    # generate sections
    sec1 = gen_sec1(dim_codes_y0, dim_codes_y1, dim_codes_cb, dim_codes_cr)
    sec2, sec2_linepos = gen_sec2(comp_codes)
    sec3, sec3_linepos = gen_sec3(comp_nums)
    sec0 = gen_sec0(sec2_linepos, sec3_linepos)
    
    # generate header
    sec1_start = 10 + len(sec0)
    sec2_start = sec1_start + len(sec1)
    sec3_start = sec2_start + len(sec2)
    if sec3_start % 4:
        padding = 4 - (sec3_start % 4)
        sec2 += bytearray(padding)
        sec3_start += padding
    frame_header = bytearray(10)
    frame_header[0] = 0b0 | level
    frame_header[1] = (level << 4) | level
    frame_header[2] = width
    frame_header[3] = height
    frame_header[4] = sec1_start >> 8
    frame_header[5] = sec1_start & 0xff
    frame_header[6] = sec2_start >> 8
    frame_header[7] = sec2_start & 0xff
    frame_header[8] = sec3_start >> 8
    frame_header[9] = sec3_start & 0xff
    
    res = frame_header+sec0+sec1+sec2+sec3
    if len(res) % 0x800:
        padding = 0x800 - (len(res) % 0x800)
        res += bytearray(padding)
        
    return res