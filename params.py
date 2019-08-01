#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

"""
Various parameters controlling SKReact, kept here for cleanliness
"""

import pandas as pd

# Reactor info filepath
REACT_FILE_PATH = "./react_p/DB2015-17.xls"

# SKReact Parameters
WIN_X = 900
WIN_Y = 1000

# Super-Kamiokande Info
SK_LAT = 36.4267    # deg
SK_LONG = 137.3104   # deg
SK_ALT = 0.370      # km
SK_FID_MASS = 22.5  # kt

# E Spectrum Hyper-parameters
E_MIN = 0 # MeV
E_MAX = 10 # MeV
E_BINS = 100
E_INTERVAL = (E_MAX-E_MIN)/E_BINS

# Reactor nu production
NU_PER_FISS = 6
NU_PER_MW = 2e17 #/s
# Using first values from PHYSICAL REVIEW D 91, 065002 (2015)
# TODO: Which values are best? Also, MOX is different for different reactors
core_types = ["PWR","BWR","LWGR","GCR","PHWR","MOX"]
__fuel_makeup_data = {"U_235": [0.538,0.538,0.538,0.538,0.560,0.560],
                "Pu_239":[0.328,0.328,0.328,0.328,0.300,0.300],
                "U_238": [0.078,0.078,0.078,0.078,0.080,0.080],
                "Pu_241":[0.056,0.056,0.056,0.056,0.060,0.060]}
FUEL_MAKEUP = pd.DataFrame(__fuel_makeup_data, index=core_types)

# 5th order polynomial coeffs and error for E_spectra 1 thru 6
# TODO: Put into dataframe
E_SPEC_N_ORDER = 5
U_235_A = [3.217, -3.111, 1.395, -3.690e-1, 4.445e-2, -2.053e-3]
U_235_DA = [4.09e-2, 2.34e-2, 4.88e-3, 6.08e-4, 7.77e-5, 6.79e-6]

PU_239_A = [6.413, -7.432, 3.535, -8.820e-1, 1.025e-1, -4.550e-3]
PU_239_DA = [4.57e-2, 2.85e-2, 6.44e-3, 9.11e-4, 1.38e-4, 1.29e-5]

U_238_A = [4.833e-1, 1.927e-1, -1.283e-1, -6.762e-3, 2.233e-3, -1.536e-4]
U_238_DA = [1.24e-1, 5.86e-2, 1.11e-2, 1.92e-3, 2.84e-4, 2.86e-5]

PU_241_A = [3.251, -3.204, 1.428, -3.675e-1, 4.254e-2, -1.896e-3]
PU_241_DA = [4.37e-2, 2.60e-2, 5.66e-3, 7.49e-4, 1.02e-4, 9.03e-6]

# Misc.
EARTH_R = 6371 # km
EARTH_R_POLAR = 6356 # km
EARTH_R_EQUATOR = 6378 # km 
