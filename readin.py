#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:06:13 2020

@author: fabiopacucci
"""

import pylab
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import stats
import random
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u



hdul = fits.open('DR7.fits')
data = hdul[1].data  # First extension is the table!
ra = data['RA']
dec = data['DEC']
z = data['REDSHIFT']
mass = data['LOGBH']

n_zbin = 50
z_min = 1.0
z_max = 5.0
z_bin = np.linspace(z_min, z_max, n_zbin)

for i in range(0, n_zbin-1):

    z_min = z_bin[i]
    z_max = z_bin[i+1]
    

    z_new = z[(z > z_min) & (z < z_max)]
    ra_new = ra[(z > z_min) & (z < z_max)]
    dec_new = dec[(z > z_min) & (z < z_max)]
    mass_new = mass[(z > z_min) & (z < z_max)]
