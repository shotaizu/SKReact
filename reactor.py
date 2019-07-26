#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

import matplotlib as plt
import pandas as pd
import os


# Getting all reactor power information into pd df
react_dat_path = "./react_p/"
react_dat = []

for r, d, f in os.walk(react_dat_path):
    for file in f:
        if ".xls" in file:
            # TODO: Handle years (rather than arbitrary index)
            react_dat.append(pd.read_excel(react_dat_path + file,
                index_col=1, 
                names=["Country Code",
                    "Reactor Core",
                    "Latitude",
                    "Longitude",
                    "Core Type",
                    "Use MOX?",
                    "Thermal Power",
                    "LF_Jan",
                    "LF_Feb",
                    "LF_Mar",
                    "LF_Apr",
                    "LF_May",
                    "LF_Jun",
                    "LF_Jul",
                    "LF_Aug",
                    "LF_Sep",
                    "LF_Oct",
                    "LF_Nov",
                    "LF_Dec"]))

for dat in react_dat:
    print(dat.loc["TAKAHAMA-3"])
