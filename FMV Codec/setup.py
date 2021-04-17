# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 08:30:57 2020

@author: nanashi
"""

from distutils.core import setup

setup (name = 'ss_grandia_codec',
       version = '0.6',
       description = 'Sega Saturn Grandia FMV codec',
       packages=['ss_grandia_codec'],
       package_data={'ss_grandia_codec': ['order.npy', 'quantize.npy']},
       install_requires='numpy >=1.19.3',
       python_requires='>=3.8')