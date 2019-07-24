#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

import pandas as pd
import os


# Getting all reactor power information into pd df
react_dat_path = "./react_p/"
react_dat = []

for r, d, f in os.walk(react_dat_path):
    for file in f:
        if ".xls" in file:
            # TODO: Handle years (rather than arbitrary index)
            react_dat.append(pd.read_excel(react_dat_path + file,index_col=1))

for dat in react_dat:
    print(dat.loc["TAKAHAMA-3",:])
