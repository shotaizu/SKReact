#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
    Produces a smearing matrix based on WIT energy reconstruction resolution
    for the SKReact app
"""

from params import *
import pandas as pd
import numpy as np
import math

def gaussian(x, mu, sig, c=1):
    return c*np.exp(-(x - mu)**2) / (2*(sig**2))

class Smear:
    def __init__(self, filename):
        wit_dat = pd.read_csv(filename, index_col = "e")
        # Create series of energies
        energies_df = pd.DataFrame(np.nan, 
            index = ENERGIES, 
            columns = ["c","mu","sig","eff"])
        energies_df = energies_df[~energies_df.index.isin(wit_dat.index)]
        # Concat with wit dat
        full_dat = pd.concat([wit_dat,energies_df])
        full_dat.sort_index(inplace=True)
        # Interpolate to fill the new values between and AFTER WIT points
        full_dat.interpolate(method = "linear",
            limit_direction = "forward",
            inplace = True)
        # Check the shape of the matrix works
        if(full_dat.shape[0] != E_BINS):
            print("SHAPE ISSUE IN SMEARING MATRIX")
            print("Check for floating point errors in smearing data.")

        # Calc gaussian for each row with SMEAR_BINS bins
        # Modify area according to efficiency
        gauss_list = []
        for row in full_dat.itertuples():
            # Check if it is below the WIT smear data
            # Assume it won't be detected at all if so
            if(math.isnan(row.mu)):
                smear_gauss = [0] * SMEAR_BINS
            else:
                smear_gauss = [
                    gaussian(energy, 
                        row.mu,
                        row.sig,
                        row.eff*row.c)
                    for energy in SMEAR_ENERGIES]
                
            gauss_list.append(smear_gauss)

        self.smear_mat = np.vstack(gauss_list)
        self.inverse_smear = np.linalg.inv(self.smear_mat)
        return

    """
        Takes a pandas Series input spectrum and multiplies by the smearing 
        matrix to produce a smeared spectrum of same length
    """
    def smear(self, spec):
       return(np.matmul(spec.to_numpy(),self.smear_mat))

    """
        Calculates inverse smearing matrix for unfolding from the already
        calculated matrix
    """
    def inverse_smear(self, spec):
       return(np.matmul(spec.to_numpy(),self.inverse_mat))