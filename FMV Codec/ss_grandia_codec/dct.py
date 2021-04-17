# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 14:29:45 2021

@author: nanashi
"""

import numpy as np
    
dct_mat_1_fix = np.array(
    [0x5a82, 0x7d8a, 0x7642, 0x6a6e, 0x5a82, 0x471d, 0x30fc, 0x18f9,
     0x5a82, 0x6a6e, 0x30fc,-0x18f9,-0x5a82,-0x7d8a,-0x7642,-0x471d,         
     0x5a82, 0x471d,-0x30fc,-0x7d8a,-0x5a82, 0x18f9, 0x7642, 0x6a6e, 
     0x5a82, 0x18f9,-0x7642,-0x471d, 0x5a82, 0x6a6e,-0x30fc,-0x7d8a, 
     0x5a82,-0x18f9,-0x7642, 0x471d, 0x5a82,-0x6a6e,-0x30fc, 0x7d8a, 
     0x5a82,-0x471d,-0x30fc, 0x7d8a,-0x5a82,-0x18f9, 0x7642,-0x6a6e, 
     0x5a82,-0x6a6e, 0x30fc, 0x18f9,-0x5a82, 0x7d8a,-0x7642, 0x471d, 
     0x5a82,-0x7d8a, 0x7642,-0x6a6e, 0x5a82,-0x471d, 0x30fc,-0x18f9
     ], dtype="int64").reshape((8,8))

dct_mat_21_fix = np.array(
    [[0x5a82, 0x5a82, 0x5a82, 0x5a82, 0x5a82, 0x5a82, 0x5a82, 0x5a82]], dtype='int64')
dct_mat_22_fix = np.array(
    [[0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200],
     [0x7d8a80, 0x6a6080, 0x471d00, 0x192200,-0x192200,-0x471d00,-0x6a6080,-0x7d8a80]], dtype='int64')
dct_mat_23_fix = np.array(
    [[0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200],
     [0x6a6d80, 0x7d8a80, 0x192200, 0x471d00,-0x471d00,-0x192200,-0x7d8a80,-0x6a6d80],
     [0x30fc00, 0x764200,-0x764200,-0x30fc00,-0x30fc00,-0x764200, 0x764200, 0x30fc00]], dtype='int64')
dct_mat_24_fix = np.array(
    [[ 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200, 0x5a8200],
     [ 0x6a6d80, 0x7d8a80, 0x192200, 0x471d00,-0x471d00,-0x192200,-0x7d8a80,-0x6a6d80],
     [ 0x30fc00, 0x764200,-0x764200,-0x30fc00,-0x30fc00,-0x764200, 0x764200, 0x30fc00],
     [-0x18f900, 0x6a6e00,-0x471d00, 0x7d8a00,-0x7d8a00, 0x471d00,-0x6a6e00,-0x18f900]], dtype='int64')
    

def _gen_dct_mat():
    res = np.zeros((8,8))
    for j in range(0,8):
        for k in range(0,8):
            if j == 0:
                en = np.sqrt(2)/2
            else:
                en = 1
            res[j,k] = en*np.cos((j*(2*k+1)*np.pi)/16)
    return 0.5*res

dct_mat1 = _gen_dct_mat()
dct_mat2 = _gen_dct_mat().T

def _idct_21(block):
    
    res = np.dot(block,dct_mat_21_fix)
    return np.right_shift(res, 32)

def _idct_22(block):
    
    res1 = np.dot(block[:,0:1], dct_mat_22_fix[0:1,:])
    res2 = np.dot(block[:,1:2], dct_mat_22_fix[1:2,:])
    res = np.right_shift(res1, 32) + np.right_shift(res2, 32)
    res = np.right_shift(res, 8)
    return res

def _idct_23(block):
    
    res1 = np.dot(block[:,0:1], dct_mat_23_fix[0:1,:])
    res2 = np.dot(block[:,1:2], dct_mat_23_fix[1:2,:])
    res3 = np.dot(block[:,2:3], dct_mat_23_fix[2:3,:])
    res = np.right_shift(res1, 32) + np.right_shift(res2, 32) + np.right_shift(res3, 32)
    res = np.right_shift(res, 8)
    return res

def _idct_24(block):
    
    res1 = np.dot(block[:,0:1], dct_mat_24_fix, block[0:1,:])
    res2 = np.dot(block[:,1:2], dct_mat_24_fix, block[1:2,:])
    res3 = np.dot(block[:,2:3], dct_mat_24_fix, block[2:3,:])
    res4 = np.dot(block[:,3:4], dct_mat_24_fix, block[3:4,:])
    res = np.right_shift(res1, 32) + np.right_shift(res2, 32) + \
        np.right_shift(res3, 32) + np.right_shift(res4, 32)
    res = np.right_shift(res, 8)
    return res

def _idct_25(block):
    
    res = np.empty((8,8), dtype="int64")
    for i in range(0,8):
        fp = 0x5a8200*(block[i,0] + block[i,4]) >> 32
        fy = 0x5a8200*(block[i,0] - block[i,4]) >> 32
        
        wtf = fy +((0x471d00*block[i,1])>>32) -((0x30fc00*block[i,2])>>32)                  # The fuck is this?!
        wtf = (0x7d8a00*wtf)>>32
        
        res[i,0] = fy +((0x6a6d80*block[i,1])>>32) +((0x30fc00*block[i,2])>>32) -((0x18F900*block[i,3])>>32)
        res[i,1] = fp +((0x7d8a00*block[i,1])>>32) +((0x764200*block[i,2])>>32) +((0x6a6e00*block[i,3])>>32)
        res[i,2] = fp +((0x192200*block[i,1])>>32) -((0x764200*block[i,2])>>32) -((0x471d00*block[i,3])>>32)
        res[i,3] = fy -((0x30fc00*block[i,1])>>32) +wtf
        res[i,4] = fy -((0x30fc00*block[i,1])>>32) -wtf
        res[i,5] = fp -((0x192200*block[i,1])>>32) -((0x764200*block[i,2])>>32) +((0x471d00*block[i,3])>>32)
        res[i,6] = fp -((0x7d8a00*block[i,1])>>32) +((0x764200*block[i,2])>>32) -((0x6a6e00*block[i,3])>>32)
        res[i,7] = fy -((0x6a6d80*block[i,1])>>32) +((0x30fc00*block[i,2])>>32) +((0x18F900*block[i,3])>>32)
        
    return np.right_shift(res, 8)

def _idct_25c(block):
    
    res = np.empty((8,8), dtype="int64")
    for i in range(0,8):
        fp = 0x5a8200*(block[i,0] + block[i,4]) >> 32
        fy = 0x5a8200*(block[i,0] - block[i,4]) >> 32
        
        res[i,0] = fp +((0x7d8a00*block[i,1])>>32) +((0x764200*block[i,2])>>32) +((0x6a6e00*block[i,3])>>32)
        res[i,1] = fy +((0x6a6d80*block[i,1])>>32) +((0x30fc00*block[i,2])>>32) -((0x18F900*block[i,3])>>32)
        res[i,2] = fy +((0x471d00*block[i,1])>>32) -((0x30fc00*block[i,2])>>32) +((0x7d8a00*block[i,3])>>32)
        res[i,3] = fp +((0x192200*block[i,1])>>32) -((0x764200*block[i,2])>>32) -((0x471d00*block[i,3])>>32)
        res[i,4] = fp -((0x192200*block[i,1])>>32) -((0x764200*block[i,2])>>32) +((0x471d00*block[i,3])>>32)
        res[i,5] = fy -((0x471d00*block[i,1])>>32) -((0x30fc00*block[i,2])>>32) -((0x7d8a00*block[i,3])>>32)
        res[i,6] = fy -((0x6a6d80*block[i,1])>>32) +((0x30fc00*block[i,2])>>32) +((0x18F900*block[i,3])>>32)
        res[i,7] = fp -((0x7d8a00*block[i,1])>>32) +((0x764200*block[i,2])>>32) -((0x6a6e00*block[i,3])>>32)
        
    return np.right_shift(res, 8)

def _idct_28(block):
    
    res = np.empty((8,8), dtype="int64")
    for i in range(0,8):
        fp = 0x5a8200*(block[i,0] + block[i,4]) >> 32
        fy = 0x5a8200*(block[i,0] - block[i,4]) >> 32
        fb = 0x764200*block[i,2] >> 32
        fb += 0x30fc00*block[i,6] >> 32
        fgo = (0x30fc00*block[i,2]) >> 32
        fgo -= (0x764200*block[i,6]) >> 32
        fgr = 0xb50400*(block[i,3] + block[i,5]) >> 32
        fo = 0xb50400*(block[i,3] - block[i,5]) >> 32
        
        scal1 = block[i,1] >> 16
        scal1 <<= 8
        scal7 = block[i,7] >> 16
        scal7 <<= 8
        
        fgrp = scal1 + fgr
        fgrm = scal1 - fgr
        fop = scal7 + fo
        fom = scal7 - fo
        
        res[i,0] = fp +fb  +((0x7d8a8000*fgrp)>>32) +((0x19220000*fop)>>32)
        res[i,1] = fy +fgo +((0x6a6d8000*fgrm)>>32) -((0x471d0000*fom)>>32)
        res[i,2] = fy -fgo +((0x471d0000*fgrm)>>32) +((0x6a6d8000*fom)>>32)
        res[i,3] = fp -fb  +((0x19220000*fgrp)>>32) -((0x7d8a8000*fop)>>32)
        res[i,4] = fp -fb  -((0x19220000*fgrp)>>32) +((0x7d8a8000*fop)>>32)
        res[i,5] = fy -fgo -((0x471d0000*fgrm)>>32) -((0x6a6d8000*fom)>>32)
        res[i,6] = fy +fgo -((0x6a6d8000*fgrm)>>32) +((0x471d0000*fom)>>32)
        res[i,7] = fp +fb  -((0x7d8a8000*fgrp)>>32) -((0x19220000*fop)>>32)
        
    return np.right_shift(res, 8)

def dct_f(block):
    return np.dot(np.dot(dct_mat1, block), dct_mat2)

def idct_f(block):
    return np.dot(np.dot(dct_mat2, block), dct_mat1)

def idct_i(block):
    
    rows = block.shape[0]
    step_1_res = np.dot(dct_mat_1_fix[:,0:rows], block)
    
    cols = step_1_res.shape[1]
    
    if cols == 1:
        return _idct_21(step_1_res)
    elif cols == 2:
        return _idct_22(step_1_res)
    elif cols == 3:
        return _idct_23(step_1_res)
    elif cols == 4:
        return _idct_24(step_1_res)
    elif cols == 5:
        return _idct_25c(step_1_res)
    elif cols == 8:
        return _idct_28(step_1_res)
    else:
        temp = np.zeros((8,8), dtype='int64')
        temp[:,:cols] = step_1_res
        return _idct_28(temp)