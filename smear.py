#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
    Produces a smearing matrix based on WIT energy reconstruction resolution
    for the SKReact app
"""

from params import *
import pandas as pd
import numpy as np

def gaussian(x, mu, sig, c=1):
    return c*np.exp(-(x - mu)**2) / (2*(sig**2))

class Smear:
    def __init__(self, filename):
        wit_dat = pd.read_csv(filename, index_col = "e")
        for energy in ENERGIES:
            row_gauss = []
            for smear_energy in SMEAR_ENERGIES:
                # WORK OUT WHEN TO CALL NEXT ON ITERATOR
                # GO TO NEXT ROW WHEN BEYONG THE MIDPOINT BETWEEN
                # TWO ROWS IN WIT DATA
                wit_dat_row = next(wit_dat.itertuples())[1]
                # Calculate c from efficiency in file
                nearest_e = energy #CHANGE THIS TO FIND THE NEAREST BIN
                    # IN THE SMEAR FILE TO CALCULATE THE GAUSSIAN FOR 
                row_gauss.append(gaussian(smear_energy,
                    wit_dat["mu"][nearest_e],
                    wit_dat["sig"]))

    """
        Takes a pandas Series input spectrum and multiplies by the smearing 
        matrix to produce a smeared spectrum of same length
    """
    def smear(self, int_spec):
       return(np.matmul(int_spec.to_numpy(),self.smear_mat))