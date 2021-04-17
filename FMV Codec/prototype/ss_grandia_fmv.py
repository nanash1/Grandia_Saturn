# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:25:44 2020

@author: nanashi
"""

import numpy as np
import copy as cpy
from scipy.optimize import curve_fit

def gen_dct_mat1():
    
    res = np.zeros((8,8))
    
    for j in range(0,8):
        for k in range(0,8):
            if j == 0:
                en = np.sqrt(2)/2
            else:
                en = 1
                
            res[j,k] = en*np.cos((j*(2*k+1)*np.pi)/16)
            if res[j,k] < 0:
                res[j,k] += 2
            
    return (0.5*res).T

def gen_dct_mat2():
    
    res = np.zeros((8,8))
    
    for j in range(0,8):
        for k in range(0,8):
            if j == 0:
                en = np.sqrt(2)/2
            else:
                en = 1
                
            res[j,k] = en*np.cos((j*(2*k+1)*np.pi)/16)
            
    return 0.5*res

def fix_to_float(x, n):                                                         # Test function unrelated to decoding algorithm
    return x.astype(float) / (1 << n)

class section0:
    """
    Contains the information where each individually encoded line starts in
    section 2 & 3
    """
    def __init__(self, frame):
        self._line_start = []
        frame_pos = 10
        for i in range(0,11):                                                   # 11 vertical block columns per frame
            sec2_start = int().from_bytes(frame[frame_pos:frame_pos+2], "big")
            frame_pos += 2
            sec3_start = int().from_bytes(frame[frame_pos:frame_pos+2], "big")
            frame_pos += 2
            sec2_bit_offset = frame[frame_pos]
            frame_pos += 1
            sec3_bit_offset = frame[frame_pos]
            frame_pos += 1
            self._line_start.append((sec2_start, sec2_bit_offset, sec3_start, sec3_bit_offset))
            
class section1:
    """
    Contains shape codes for each macro block
    
    3 MSBs are columns
    3 LSBs are rows
    
    
    Decoded Y channel data stored @60B9000, pointer stored @602F0EC
    Decoded Cb channel data stored @60B9F20, pointer stored @602F0F0
    Decoded Cr channel data stored @60BA2E8, pointer stored @602F0F4
    """
    def __init__(self, frame):
        
        """
        Read the look up table located @6036CD0 to find the ordering of the elements
        in a macro block. Probably used for the zig zag traversal.
        """
        with open("Intro_HWRAM_dump3.bin", "rb") as ifile:
            ifile.seek(0x377F1)
            
            order_lut = []
            
            sublist_pos = ifile.tell()
            for i in range(0, 64):
                
                ifile.seek(sublist_pos)
                elem_start = int().from_bytes(ifile.read(3), "big")
                ifile.seek(1,1)
                if i == 63:
                    elem_end = 0x376F0
                else:
                    elem_end = int().from_bytes(ifile.read(3), "big")
                    sublist_pos = ifile.tell() - 3
                
                subelems = [0]*((elem_end-elem_start) >> 1)
                ifile.seek(elem_start)
                idx = 0
                while elem_start < elem_end:
                    subelems[idx] = int().from_bytes(ifile.read(2), "big") >> 2
                    idx += 1
                    elem_start += 2
                order_lut.append(subelems)
                
        self._order_lut = order_lut
        
        self._width = frame[2]                                                  # stored @602F3C8 , image width / 16
        self._heigth = frame[3]                                                 # stored @602F3CC , image height / 16
        offset = int().from_bytes(frame[4:6], "big")
        
        self._data = frame[offset:offset+self._width*self._heigth*6]
        self._pos = 0                                                           # r15 + var_44 in IDA
        self._stage = 0
        self._y_codes, self._cb_codes, self._cr_codes = self._decode()
        self._y_pos = 0
        self._cb_pos = 0
        self._cr_pos = 0
        
    def _get(self):
        """
        Fetches the next 6 bits from the frame data

        Returns
        -------
        None.

        """
        if self._stage == 0:
            res = self._data[self._pos] >> 2
        elif self._stage == 1:
            res = (self._data[self._pos] & 0x3) << 4
            self._pos += 1
            res += self._data[self._pos] >> 4
        elif self._stage == 2:
            res = (self._data[self._pos] & 0xf) << 2
            self._pos += 1
            res += self._data[self._pos] >> 6
        else:
            res = self._data[self._pos] & 0x3f
            self._pos += 1
        
        self._stage += 1
        self._stage &= 0x3
        
        return res
    
    def _align(self):
        if self._stage > 0:
            self._stage = 0
            self._pos += 2
            
    def _decode(self):
        
        y_chan = np.zeros(self._heigth*self._width*4, dtype=int)
        cb_chan = np.zeros(self._heigth*self._width, dtype=int)
        cr_chan = np.zeros(self._heigth*self._width, dtype=int)
        
        height_idx = 0
        
        while height_idx < self._heigth:
        
            width_idx = 0
            
            while width_idx < self._width:
                y_chan[88*height_idx+width_idx*2] = self._get()
                y_chan[88*height_idx+width_idx*2+1] = self._get()
                y_chan[(88*height_idx+44)+width_idx*2] = self._get()
                y_chan[(88*height_idx+44)+width_idx*2+1] = self._get()
                cb_chan[22*height_idx+width_idx] = self._get()
                cr_chan[22*height_idx+width_idx] = self._get()
                
                self._align()
                width_idx += 1
                
            height_idx += 1
            
        return y_chan, cb_chan, cr_chan
    
    def _get_elem_num(self, code):                                              # game actually uses a lut that starts @6036C90
        return ((code & 0x7) + 1) * ((code >> 3) + 1)
    
    def get_dim(self, code):
        """
        Returns dimension of the sub macro block. The result is swapped
        compared to numpy's notation because the game fills columns first!

        Parameters
        ----------
        code : int
            dimension code from section 1.

        Returns
        -------
        List
            (columns, rows).
        """                                                   
        return ((code >> 3) + 1, (code & 0x7) + 1)
    
    def get_y(self):
        code = self._y_codes[self._y_pos]
        num = self._get_elem_num(code)
        order = cpy.copy(self._order_lut[code])
        self._y_pos += 1
        return code, num, order
    
    def get_cb(self):
        code = self._y_codes[self._cb_pos]
        num = self._get_elem_num(code)
        order = cpy.copy(self._order_lut[code])
        self._cb_pos += 1
        return code, num, order
    
    def get_cr(self):
        code = self._y_codes[self._cr_pos]
        num = self._get_elem_num(code)
        order = cpy.copy(self._order_lut[code])
        self._cr_pos += 1
        return code, num, order
    
    def reset(self):
        self._y_pos = 0
        self._cb_pos = 0
        self._cr_pos = 0
    
class section2:
    """
    Contains the instruction codes that tell us how to decode the macro block
    elements.
    """
    def __init__(self, frame, sec2_offset):
        
        lut_repeats = [0, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01]
        lut_seed = [0]*32                                                           # simulated memory @6031454
        lut_table = [0]*256                                                         # simulated memory @6031480
    
        idx = 0
        sec2_pos = sec2_offset
        while idx < 32:        
            lut_seed[idx] = (frame[sec2_pos] >> 4) & 0xf
            idx += 1
            
            lut_seed[idx] = frame[sec2_pos] & 0xf
            idx += 1
            
            sec2_pos += 1
            
        idx = 0
        while idx < 256:
            
            max_ind = lut_seed.index(max(lut_seed))
            max_elem = lut_seed[max_ind]
            lut_seed[max_ind] = 0
            new_table_elem = (max_elem << 4) | (max_ind << 9)
            
            repeats = lut_repeats[max_elem]
            
            while repeats:
                lut_table[idx] = new_table_elem
                repeats -= 1
                idx += 1
        self._frame = frame
        self._lut = lut_table
        self._pos = sec2_offset + 16
        self._bit_pos = 0
        self._window = self._frame[self._pos]
        self._pos += 1
        
    def get(self):
        
        res = self._lut[self._window]
        shift = (res >> 4) & 0xf
        res >>= 9
        
        if self._bit_pos + shift > 7:
            self._window <<= 8 - self._bit_pos
            self._window &= (0xff << (8 - self._bit_pos)) & 0xff
            self._window |= self._frame[self._pos] & (2**(8-self._bit_pos) - 1)
            self._pos += 1
            self._bit_pos += shift
            self._bit_pos &= 0x7
            self._window <<= self._bit_pos
        else:
            self._window <<= shift
            self._bit_pos += shift
            
        self._window &= (0xff << self._bit_pos) & 0xff
        self._window |= self._frame[self._pos] >> (8 - self._bit_pos)
        
        return res
        
class section3:
    """
    Contains the the acutal macro block coefficient data. However, when the
    coefficients are simple multiples of the quantization factors, these numbers
    aren't used.
    """
    def __init__(self, frame, sec3_offset):
        self._frame = frame
        self._pos = sec3_offset
        self._bit_pos = 8
        self._window = 0
        
    def get(self, bits):
        
        res = 0
        
        while bits:
        
            new_bits = self._frame[self._pos] & (2**self._bit_pos - 1)
            
            if bits > self._bit_pos:
                res <<= self._bit_pos
                bits -= self._bit_pos
                self._bit_pos = 8
                self._pos += 1
            else:
                res <<= bits
                new_bits >>= self._bit_pos - bits
                self._bit_pos -= bits
                bits = 0
                
            res |= new_bits
            
        return res
    
class instruction_decoder:
    
    """
    The instruction table decodes the instruction code. Starts @6039AF4. 
    
    Each Element consists of these instructions
        - number of zeros to write at the begining
        - how many bits of section 3 to read
        - offset to add to the section 3 bits
        - quantization factor to multiply look up table data with
        
    table structure:
        (number of zeros, sec3_bits, sec3_offset, quant_fac)
    """
    _instruction_table =  [(0, 0,  0,  1),                              #0
                           (1, 0,  0,  1),                              #1
                           (2, 0,  0,  1),                              #2
                           (3, 0,  0,  1),                              #3
                           (4, 0,  0,  1),                              #4
                           (0, 0,  0, -1),                              #5
                           (1, 0,  0, -1),                              #6
                           (2, 0,  0, -1),                              #7
                           (3, 0,  0, -1),                              #8
                           (4, 0,  0, -1),                              #9
                           (0, 0,  0,  2),                              #10
                           (0, 0,  0, -2),                              #11
                           (1, 0,  0,  2),                              #12
                           (1, 0,  0, -2),                              #13
                           (0, 0,  0,  3),                              #14
                           (0, 0,  0, -3),                              #15
                           (0, 0,  0,  4),                              #16
                           (0, 0,  0, -4),                              #17
                           (0, 1,  5,  1),                              #18
                           (0, 1,  5, -1),                              #19
                           (0, 2,  7,  1),                              #20
                           (0, 2,  7, -1),                              #21
                           (0, 3, 11,  1),                              #22
                           (0, 3, 11, -1),                              #23
                           (0, 4, 19,  1),                              #24
                           (0, 4, 19, -1),                              #25
                           (0, 8, 35,  1),                              #26
                           (0, 8, 35, -1)]                              #27
    """
    Used to generate macro block coefficients from instruction codes
    """
    def __init__(self, sec3_data):
        
        self._sec3_data = sec3_data
        """
        This look up table contains the quantization factors
        """
        with open("Intro_HWRAM_dump3.bin", "rb") as ifile:
            ifile.seek(0x376F1)
            
            self._quant_facs = []
            
            sublist_pos = ifile.tell()
            for i in range(0, 64):
                ifile.seek(sublist_pos)
                elem_start = int().from_bytes(ifile.read(3), "big")
                ifile.seek(1,1)
                elem_end = int().from_bytes(ifile.read(3), "big")
                sublist_pos = ifile.tell() - 3
                
                elem_num = ((i & 7) + 1) * ((i>>3) + 1)
                subelems = []
                ifile.seek(elem_start)
                idx = 0
                while elem_start < elem_end:
                    
                    subsubelems = []
                    for i in range(0, elem_num):
                        subsubelems.append(int().from_bytes(ifile.read(1), "big"))
                        idx += 1
                        elem_start += 1
                        
                    subelems.append(subsubelems)
                self._quant_facs.append(subelems)
                
        self._quant_sel = None
        self._quant_pos = 0
        self._elem_cntr = 0
        self._elem_num = 0
        
    def set_start(self, sec1_code, quant_sel):
        self._quant_sel = self._quant_facs[sec1_code][quant_sel]
        self._quant_pos = 0
        self._elem_cntr = 0
                
    def decode_instr(self, instruction_code): 
        
        rtn = []
        
        if instruction_code == 30:                                              # send "end of block" instruction
            rtn.append("EOB")
            
        elif instruction_code == 29:                                            # interpret section 3 data as number of 0s to write
            sec3_elem = self._sec3_data.get(3)
            sec3_elem += 2
            
            while sec3_elem:
                rtn.append(0)
                self._quant_pos += 1
                self._elem_cntr += 1
                sec3_elem -= 1
                
        elif instruction_code == 28:                                            # write only one zero and don't use lut data
            rtn.append(0)
            self._quant_pos += 1
            self._elem_cntr += 1
            
        else:
            instruction = self._instruction_table[instruction_code]             # use instruction table to decode instruction code
            
            for i in range(0, instruction[0]):                                  # write zeros
                rtn.append(0)
                self._quant_pos += 1
                self._elem_cntr += 1
                
            decoded_elem = 1
            if instruction[1] > 0:                                              # read section 3 bits
                if instruction[1] < 5:
                    decoded_elem = self._sec3_data.get(instruction[1])
                elif self._elem_cntr == 0:
                    decoded_elem = self._sec3_data.get(11)
                elif self._elem_cntr < 3:
                    decoded_elem = self._sec3_data.get(10)
                else:
                    decoded_elem = self._sec3_data.get(8)
                    
                decoded_elem += instruction[2]
                
            quant_fac = self._quant_sel[self._quant_pos]                        # get quantization factor
            self._quant_pos += 1
            decoded_elem *= quant_fac * instruction[3]
            
            rtn.append(decoded_elem)                                            # write macro block element
            self._elem_cntr += 1
            
        return rtn
    
class dct:
    
    '''
    DCT matrix where all negative entries are offset by +1
    '''
    dct_mat_1_fix_ref = np.array(
        [0x5a82, 0x7d8a, 0x7642, 0x6a6e, 0x5a82, 0x471d, 0x30fc, 0x18f9,    # located @6039844
         0x5a82, 0x6a6e, 0x30fc, 0xe707, 0xa57e, 0x8276, 0x89be, 0xb8e3,         
         0x5a82, 0x471d, 0xcf04, 0x8276, 0xa57e, 0x18f9, 0x7642, 0x6a6e, 
         0x5a82, 0x18f9, 0x89be, 0xb8e3, 0x5a82, 0x6a6e, 0xcf04, 0x8276, 
         0x5a82, 0xe707, 0x89be, 0x471d, 0x5a82, 0x9592, 0xcf04, 0x7d8a, 
         0x5a82, 0xb8e3, 0xcf04, 0x7d8a, 0xa57e, 0xe707, 0x7642, 0x9592, 
         0x5a82, 0x9592, 0x30fc, 0x18f9, 0xa57e, 0x7d8a, 0x89be, 0x471d, 
         0x5a82, 0x8276, 0x7642, 0x9592, 0x5a82, 0xb8e3, 0x30fc, 0xe707
         ], dtype="int64").reshape((8,8))
    
    dct_mat_1_fix = np.array(
        [0x5a82, 0x7d8a, 0x7642, 0x6a6e, 0x5a82, 0x471d, 0x30fc, 0x18f9,    # located @6039844
         0x5a82, 0x6a6e, 0x30fc,-0x18f9,-0x5a82,-0x7d8a,-0x7642,-0x471d,         
         0x5a82, 0x471d,-0x30fc,-0x7d8a,-0x5a82, 0x18f9, 0x7642, 0x6a6e, 
         0x5a82, 0x18f9,-0x7642,-0x471d, 0x5a82, 0x6a6e,-0x30fc,-0x7d8a, 
         0x5a82,-0x18f9,-0x7642, 0x471d, 0x5a82,-0x6a6e,-0x30fc, 0x7d8a, 
         0x5a82,-0x471d,-0x30fc, 0x7d8a,-0x5a82,-0x18f9, 0x7642,-0x6a6e, 
         0x5a82,-0x6a6e, 0x30fc, 0x18f9,-0x5a82, 0x7d8a,-0x7642, 0x471d, 
         0x5a82,-0x7d8a, 0x7642,-0x6a6e, 0x5a82,-0x471d, 0x30fc,-0x18f9
         ], dtype="int64").reshape((8,8))
    
    '''
    Regular DCT matrix being used (just for reference)
    '''
    dct_mat_2_fix = np.array(
        [0x5a820000, 0x5a820000, 0x5a820000, 0x5a820000, 0x5a820000, 0x5a820000, 0x5a820000, 0x5a820000, 
         0x7d8a8000, 0x6a6d8000, 0x471d0000, 0x19220000,-0x19220000,-0x471d0000,-0x6a6d8000,-0x7d8a8000,
         0x76420000, 0x30fc0000,-0x30fc0000,-0x76420000,-0x76420000,-0x30fc0000, 0x30fc0000, 0x76420000,
         0x6a8a5132,-0x18f877c2,-0x7d89b2aa,-0x46ff7422, 0x46ff7422, 0x7d89b2aa, 0x18f877c2,-0x6a8a5132,
         0x5a820000,-0x5a820000,-0x5a820000, 0x5a820000, 0x5a820000,-0x5a820000,-0x5a820000, 0x5a820000,
         0x46ff7422,-0x7d89b2aa, 0x18f877c2, 0x6a8a5132,-0x6a8a5132,-0x18f877c2, 0x7d89b2aa,-0x46ff7422,
         0x30fc0000,-0x76420000, 0x76420000,-0x30fc0000,-0x30fc0000, 0x76420000,-0x76420000, 0x30fc0000,
         0x19220000,-0x471d0000, 0x6a6d8000,-0x7d8a8000, 0x7d8a8000,-0x6a6d8000, 0x471d0000,-0x19220000
         ], dtype="int64").reshape((8,8))
    
    dct_2_2 = np.array(
        [0x7D8A80, 0x6A6080, 0x471D00, 0x192200,-0x192200,-0x471D00,-0x6A6080,-0x7D8A80], dtype='int64')
    
    dct_2_stage_18 = np.array(
        [[ 1, 0x764200, 0x30fc00, 1, 1],
         [-1, 0x30fc00,-0x764200,-1,-1],
         [-1, 0x30fc00,-0x764200,-1,-1],
         [ 1, 0x764200, 0x30fc00, 1, 1],
         [ 1, 0x764200, 0x30fc00, 1, 1],
         [-1, 0x30fc00,-0x764200,-1,-1],
         [-1, 0x30fc00,-0x764200,-1,-1],
         [ 1, 0x764200, 0x30fc00, 1, 1]], dtype="int64")
       
    dct_2_stage_28 = np.array(
        [[1, 1, 0x7d8a8000, 0x19220000],
         [1, 1, 0x6a6d8000,-0x471d0000],
         [1,-1, 0x471d0000, 0x6a6d8000],
         [1,-1, 0x19220000,-0x7d8a8000],
         [1,-1,-0x19220000, 0x7d8a8000],
         [1,-1,-0x471d0000,-0x6a6d8000],
         [1, 1,-0x6a6d8000, 0x471d0000],
         [1, 1,-0x7d8a8000,-0x19220000]], dtype="int64")
    
    def idct(self, block):
        
        step_1_res = self._idct_1(block)
        
        dim = step_1_res.shape[1]
        
        if dim == 1:
            return self._idct_21(step_1_res)
        elif dim == 2:
            return self._idct_22(step_1_res)
        else:
            return self._idct_28(step_1_res)
    
    def _idct_1(self, block):
        
        rows = block.shape[0]
        
        return np.dot(self.dct_mat_1_fix[:,0:rows], block)
    
    def _idct_21(self, block):
        
        decoded_block = np.empty((8,8), dtype="int64")
        for i in range(0,8):
            decoded_block[:,i] = 0x5a82 * block.reshape(8)
        return np.right_shift(decoded_block, 32)
    
    def _idct_22(self, block):
        
        decoded_block = np.empty((8,8), dtype="int64")
        facs = np.zeros((8,2), dtype='int64')
        for i in range(0,8):
            facs[:,0] = 0x5a8200*block[:,0]
            facs[:,1] = self.dct_2_2[i]*block[:,1]
            facs[:,0] >>= 32
            facs[:,1] >>= 32
            
            decoded_block[:,i] = facs.sum(axis=1)
            decoded_block[:,i] >>= 8
            
        return decoded_block
        
    
    def _idct_28(self, block):
        
        decoded_block = np.empty((8,8), dtype="int64")
        facs = np.zeros((8,4), dtype='int64')
        
        '''
        prescale inputs
        '''
        presc = block       
        presc[:,1] >>= 16
        presc[:,1] <<= 8
        presc[:,7] >>= 16
        presc[:,7] <<= 8
        
        '''
        calculate 2nd iDCT step
        '''
        for j in range(0,8):
            
            facs[:,0] = 0x5a8200*(presc[:,0]+self.dct_2_stage_18[j,0]*presc[:,4]) >> 32
            facs[:,1] = ((presc[:,2]*self.dct_2_stage_18[j,1]) >> 32) + \
                ((presc[:,6]*self.dct_2_stage_18[j,2]) >> 32)
            facs[:,2] = presc[:,1]+self.dct_2_stage_18[j,3]*(0xb50400*(presc[:,3]+presc[:,5]) >> 32)
            facs[:,3] = presc[:,7]+self.dct_2_stage_18[j,3]*(0xb50400*(presc[:,3]-presc[:,5]) >> 32)
            
            for k in range(0,4):
                facs[:,k] *= self.dct_2_stage_28[j,k] 
            facs[:,2] >>= 32
            facs[:,3] >>= 32
            
            decoded_block[:,j] = facs.sum(axis=1)
            decoded_block[:,j] >>= 8
            
        return decoded_block
    
def YCbCr2RGB(y0, y1, cr, cb):
    
    def mask(x):
        if x > 0:
            return x & 0xff
        else:
            return x & 0x00
        
    """
    y stored in M1
    cr stored in M2
    cb stored in M3
    """
    y0 = int().from_bytes(y0[0:2], 'big', signed=True)
    y1 = int().from_bytes(y0[2:4], 'big', signed=True)
    y2 = int().from_bytes(y1[0:2], 'big', signed=True)
    y3 = int().from_bytes(y1[2:4], 'big', signed=True)
    cr = int().from_bytes(cr, 'big', signed=True)
    cb = int().from_bytes(cb, 'big', signed=True)
    
    y0 += 124
    y1 += 124
    y2 += 124
    y3 += 124
    
    green0 = mask(y0 + cb + cr)    
    green1 = mask(y1 + cb + cr)
    green2 = mask(y2 + cb + cr)
    green3 = mask(y3 + cb + cr)
    red0 = mask(-2*cr + y0)
    red1 = mask(-2*cr + y2) 
    blue0 = mask(-2*cb + y0)
    blue1 = mask(-2*cb + y2)
        
    rgb0 = blue0 << 16 | green0 << 8 | red0
    rgb1 = blue0 << 16 | green1 << 8 | red0
    rgb2 = blue1 << 16 | green2 << 8 | red1
    rgb3 = blue1 << 16 | green3 << 8 | red1
    
    return rgb0, rgb1, rgb2, rgb3
        

def decode_differential(macro_blocks):
    elem0 = 0
    for macro_block in macro_blocks:
        elem0 += macro_block[0,0]
        macro_block[0,0] = elem0
    
def gen_block(frame, sec1_data, sec2_data, instr_decoder, quant_sel):
    
    sec1_code, elem_num, order = sec1_data.get_y()
    macro_block = np.empty(elem_num, dtype="int64")
    instr_decoder.set_start(sec1_code, quant_sel)
    
    while True:
        instruction_code = sec2_data.get()
        block_elems = instr_decoder.decode_instr(instruction_code)
        
        for elem in block_elems:
            if elem == "EOB":
                while elem_num:
                    elem_num -= 1
                    macro_block[order.pop(0)] = 0
                return macro_block.reshape(sec1_data.get_dim(sec1_code)).T
            macro_block[order.pop(0)] = elem
            elem_num -= 1
            
def decode_macroblocks(frame, sec1_data, sec2_data, instr_decoder, quant_sel):
    
    macro_blocks = []                                                           #first block @60C56B0, second block @60C57B0, third block @☺60C58B0  
    cntr = 44    
    while cntr:
        if cntr == 5:
            test = 1
        macro_blocks.append(gen_block(frame, sec1_data, sec2_data, instr_decoder, quant_sel))
        cntr -= 1   
        
    decode_differential(macro_blocks)
    sec1_data.reset()
    
    y_blocks = []                                                               # stored @603FFC0, offset 0x2c0 to new line because 2 bytes * 22 blocks * 16 pixel
    odct = dct()                                                                # @607C7C0
    
    for i in range(0, 44):
        
        y_blocks.append(odct.idct(macro_blocks[i]))
        np.dot(gen_dct_mat1(), np.pad(macro_blocks[6], [(0,6),(0,6)]))

        step_1_res = np.empty(64, dtype="int64")
        with open("first_full_rank_frame.bin", "rb") as ifile:
            for pos in range(0,64):
                step_1_res[pos] = int().from_bytes(ifile.read(4), "big", signed=True)
        step_1_res = step_1_res.reshape((8,8))
        odct._idct_28(step_1_res)
        
def reverse_color_conv():
    
    def fit_func(X, a, d, e, f):
        r, g, b = X
        return a*b+d*g+e*r - f
    
    r = []
    g = []
    b = []
    y = []
    cb = []
    cr = []
    with open("c:/Users/nanashi/Documents/Grandia_FMV/res1.bin", "rb") as ifile:
        for i in range(0,44):
            ifile.read(1)
            b.append(int().from_bytes(ifile.read(1), "big", signed=False))
            g.append(int().from_bytes(ifile.read(1), "big", signed=False))
            r.append(int().from_bytes(ifile.read(1), "big", signed=False))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/res2.bin", "rb") as ifile:
        for i in range(0,44):
            ifile.read(1)
            b.append(int().from_bytes(ifile.read(1), "big", signed=False))
            g.append(int().from_bytes(ifile.read(1), "big", signed=False))
            r.append(int().from_bytes(ifile.read(1), "big", signed=False))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/res3.bin", "rb") as ifile:
        for i in range(0,44):
            ifile.read(1)
            b.append(int().from_bytes(ifile.read(1), "big", signed=False))
            g.append(int().from_bytes(ifile.read(1), "big", signed=False))
            r.append(int().from_bytes(ifile.read(1), "big", signed=False))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/res4.bin", "rb") as ifile:
        for i in range(0,44):
            ifile.read(1)
            b.append(int().from_bytes(ifile.read(1), "big", signed=False))
            g.append(int().from_bytes(ifile.read(1), "big", signed=False))
            r.append(int().from_bytes(ifile.read(1), "big", signed=False))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc11.bin", "rb") as ifile:
        for i in range(0,44):
            y.append(int().from_bytes(ifile.read(2), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc12.bin", "rb") as ifile:
        for i in range(0,44):
            y.append(int().from_bytes(ifile.read(2), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc13.bin", "rb") as ifile:
        for i in range(0,44):
            y.append(int().from_bytes(ifile.read(2), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc14.bin", "rb") as ifile:
        for i in range(0,44):
            y.append(int().from_bytes(ifile.read(2), 'big', signed=True))
            
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc21.bin", "rb") as ifile:
        for i in range(0,22):
            cr.append(int().from_bytes(ifile.read(4), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc22.bin", "rb") as ifile:
        for i in range(0,22):
            cr.append(int().from_bytes(ifile.read(4), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc23.bin", "rb") as ifile:
        for i in range(0,22):
            cr.append(int().from_bytes(ifile.read(4), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc24.bin", "rb") as ifile:
        for i in range(0,22):
            cr.append(int().from_bytes(ifile.read(4), 'big', signed=True))
            
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc31.bin", "rb") as ifile:
        for i in range(0,22):
            cb.append(int().from_bytes(ifile.read(4), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc32.bin", "rb") as ifile:
        for i in range(0,22):
            cb.append(int().from_bytes(ifile.read(4), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc33.bin", "rb") as ifile:
        for i in range(0,22):
            cb.append(int().from_bytes(ifile.read(4), 'big', signed=True))
    with open("c:/Users/nanashi/Documents/Grandia_FMV/mc34.bin", "rb") as ifile:
        for i in range(0,22):
            cb.append(int().from_bytes(ifile.read(4), 'big', signed=True))
            
    r = np.array(r)
    g = np.array(g)
    b = np.array(b)
    y = np.array(y)
    cr = np.array(cr)
    cb = np.array(cb)
    
    res = curve_fit(fit_func, (r, g, b), y, (1, 1, 1, 1))
    check_res = fit_func((r,g,b), 0.25, 0.5, 0.25, 124)
    check_res2 = fit_func((r,g,b), 0.114, 0.587, 0.299, 128)
        

def decomp_frame():
    
#    YCbCr2RGB(b'\xff\xbc\xff\xbc', b'\x00\x00\x00\x0c', b'\xff\xff\xff\xf0')
#    reverse_color_conv()
    
    # frame is stored in memory @60CB000
    with open("MOV20.VID", "rb") as \
        ifile:
            
        ifile.seek(26624) # goto frame 10
        frame = ifile.read(4096)
        
    test_blk = np.load('block.npy')
    
    data = frame[0]
    var1 = data >> 4                                                            # stored @602F42C
    quant_sel_lum = data & 0xf                                                  # stored @602F3D0
    data = frame[1]
    var3 = data >> 4                                                            # stored @602F3D4
    var4 = data & 0xf                                                           # stored @602F3D8
    width = frame[2]                                                            # stored @602F3C8 , image width / 16
    height = frame[3]                                                           # stored @602F3CC , image height / 16
    sec1_offset = int().from_bytes(frame[4:6], "big")
    sec2_offset = int().from_bytes(frame[6:8], "big")                           # stored @602F414
    sec3_offset = int().from_bytes(frame[8:10], "big")                          # stored @602F418
        
    sec0_data = section0(frame)
    sec1_data = section1(frame)
    sec2_data = section2(frame, sec2_offset)
    sec3_data = section3(frame, sec3_offset)
    instr_decoder = instruction_decoder(sec3_data)
                    
    decode_macroblocks(frame, sec1_data, sec2_data, instr_decoder, quant_sel_lum)    
                
                

decomp_frame()