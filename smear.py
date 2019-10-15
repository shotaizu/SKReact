#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
    Produces a smearing matrix based on WIT energy reconstruction resolution
    for the SKReact app
"""

from params import *
import pandas as pd
import numpy as np

class Smear:
    def __init__(self, filename):
        self.wit_dat = pd.read_csv(filename, index_col = "e")

    """
        Takes a pandas Series input spectrum and multiplies by the smearing 
        matrix to produce a smeared spectrum of same length
    """
    def smear(self, int_spec):
       return(np.matmul(int_spec.to_numpy(),self.smear_mat))