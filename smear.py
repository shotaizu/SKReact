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

class Smear:
    def __init__(self):
        self.wit_dat = pd.read_csv(WIT_SMEAR_FILE, index_col = "e")

    def __init__(self,filename):
        self.wit_dat = pd.read_csv(filename, index_col = "e")