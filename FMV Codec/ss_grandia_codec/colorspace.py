# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 14:21:28 2021

@author: nanashi
"""

import numpy as np

def level_clip(channel: np.ndarray, vmin, vmax) -> np.ndarray:
    """
    Truncates values that are smaller than vmin and bigger than vmax

    Parameters
    ----------
    channel : np.ndarray
        Channel data.
    vmin : float
        Minimum channel value.
    vmax : float
        Maximum channel value.

    Returns
    -------
    np.ndarray
        Shifted array.

    """
    shape = channel.shape
    channel = channel.reshape(shape[0]*shape[1])
    idx = np.where(channel < vmin)
    channel[idx] = vmin
    idx = np.where(channel > vmax)
    channel[idx] = vmax
    return channel.reshape((shape[0], shape[1]))


def scale(channel: np.ndarray, vmin, vmax) -> np.ndarray:
    
    channel *= (vmax-vmin)/255
    return channel + vmin


def subsample420_discard(chan: np.ndarray) -> np.ndarray:
    """
    4:2:0 subsampling by discarding values 1,2,3 in a 4x4 field

    Parameters
    ----------
    chan : np.ndarray
        Channel data.

    Returns
    -------
    np.ndarray
        Subsampled channel.

    """
    return chan[::2,::2]

def subsample420_avrg(chan: np.ndarray) -> np.ndarray:
    """
    4:2:0 subsampling by averaging values 0 and 2 in a 4x4 field

    Parameters
    ----------
    chan : np.ndarray
        Channel data.

    Returns
    -------
    np.ndarray
        Subsampled channel.

    """
    chan1 = chan[::2,::2]
    chan3 = chan[1::2,::2]
    res = chan1+chan3
    return res/2

def rgb2ycbcr(r, g, b, subsample="avrg", **kwargs):
    """
    Converts RGB colorspace to YCbCr

    Parameters
    ----------
    r : np.ndarray
        Red channel.
    g : np.ndarray
        Green channel.
    b : np.ndarray
        Blue channel.
    subsample : "avrg" or "discard"
        Determines how subsampling is applied.
        Average -> Average colors vertically
        Discard -> Remove every second color point
        
    Keywords
    ----------
    rclip : tuple
        Clips red channel colors in the specified boundaries
    gclip : tuple
        Clips green channel colors in the specified boundaries
    bclip : tuple
        Clips blue channel colors in the specified boundaries
    yscale : tuple
        Scale y channel to specified boundaries
    cscale : int
        Scale cb and cr channel to specified boundaries

    Returns
    -------
    y : np.ndarray
        Luminance channel.
    cb : np.ndarray
        DESCRIPTION.
    cr : np.ndarray
        DESCRIPTION.

    """
    if "rclip" in kwargs:
        r = level_clip(r, *kwargs["rclip"])
        kwargs.pop("rclip")
    if "gclip" in kwargs:
        g = level_clip(g, *kwargs["gclip"])
        kwargs.pop("gclip")
    if "bclip" in kwargs:
        b = level_clip(b, *kwargs["bclip"])
        kwargs.pop("bclip")
    
    y = 0.25*r+0.5*g+0.25*b
    cb = (y - b)/2
    cr = (y - r)/2
    
    if "yscale" in kwargs:
        y = scale(y, *kwargs["yscale"])                                         # See Gibb's Phenomenon and https://en.wikipedia.org/wiki/YCbCr
        kwargs.pop("yscale")
    if "cscale" in kwargs:
        cb = scale(cb, 0, kwargs["cscale"])
        cr = scale(cr, 0, kwargs["cscale"])
        kwargs.pop("cscale")
    y -= 128
    
    if len(kwargs) != 0:
        raise ValueError("Unrecognized keywords")
    
    if subsample == "avrg":
        cb = subsample420_avrg(cb)
        cr = subsample420_avrg(cr)
    elif subsample == "discard":
        cb = subsample420_discard(cb)
        cr = subsample420_discard(cr)
    else:
        raise ValueError("Unknown subsample mode")
    
    return y, cb, cr

def ycbcr2rgb_420(y0, y1, y2, y3, cb, cr):
    """
    Converts YCbCr colorspace to RGB with 4:2:0 subsampling

    Parameters
    ----------
    y0 : int
        Y upper left field.
    y1 : int
        Y upper right field.
    y2 : int
        Y lower left field.
    y3 : int
        Y lower right field.
    cb : int
        Cb field.
    cr : int
        Cr field.

    Returns
    -------
    rgb0 : tuple
        Upper left field.
    rgb1 : tuple
        Upper right field.
    rgb2 : tuple
        Lower left field.
    rgb3 : tuple
        Lower right field.

    """
    y0 += 124
    y1 += 124
    y2 += 124
    y3 += 124
    
    green0 = (y0 + cb + cr) & 0xff
    green1 = (y1 + cb + cr) & 0xff
    green2 = (y2 + cb + cr) & 0xff
    green3 = (y3 + cb + cr) & 0xff
    red0 = (-2*cr + y0) & 0xff
    red1 = (-2*cr + y2) & 0xff
    blue0 = (-2*cb + y0) & 0xff
    blue1 = (-2*cb + y2) & 0xff
        
    rgb0 = (red0, green0, blue0)
    rgb1 = (red0, green1, blue0)
    rgb2 = (red1, green2, blue1)
    rgb3 = (red1, green3, blue1)

    return rgb0, rgb1, rgb2, rgb3

def ycbcr2rgb_test(y0, y1, y2, y3, cb, cr):
    """
    Converts YCbCr colorspace to RGB with 4:2:0 subsampling without adhering to
    the games overflow behavior

    Parameters
    ----------
    y0 : int
        Y upper left field.
    y1 : int
        Y upper right field.
    y2 : int
        Y lower left field.
    y3 : int
        Y lower right field.
    cb : int
        Cb field.
    cr : int
        Cr field.

    Returns
    -------
    rgb0 : tuple
        Upper left field.
    rgb1 : tuple
        Upper right field.
    rgb2 : tuple
        Lower left field.
    rgb3 : tuple
        Lower right field.

    """
    y0 += 124
    y1 += 124
    y2 += 124
    y3 += 124
    
    green0 = y0 + cb + cr
    green1 = y1 + cb + cr
    green2 = y2 + cb + cr
    green3 = y3 + cb + cr
    red0 = -2*cr + y0
    red1 = -2*cr + y2
    blue0 = -2*cb + y0
    blue1 = -2*cb + y2
        
    rgb0 = (red0, green0, blue0)
    rgb1 = (red0, green1, blue0)
    rgb2 = (red1, green2, blue1)
    rgb3 = (red1, green3, blue1)

    return rgb0, rgb1, rgb2, rgb3