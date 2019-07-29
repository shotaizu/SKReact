#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

""" 
    A GUI program to generate information regarding
    reactor neutrinos in Super-Kamiokande: production,
    oscillation, interaction and detection.
"""

from reactor import Reactor
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Getting all reactor power information into pd df
react_file_path = "./react_p/DB2015-17.xls"
react_dat = pd.read_excel(react_file_path, header=0)

# List of reactors
reactors = []
for index, data in react_dat.loc[(react_dat["Country Code"] == "JP") | (react_dat["Country Code"] == "KR")].iterrows():
    reactors.append(Reactor(
        data["Country Code"],
        data["Core Name"],
        data["Latitude"],
        data["Longitude"],
        data["Core Type"],
        data["Use Mox?"],
        data["Thermal Power"],
        data[7:]))

for reactor in reactors:
    print(reactor.n_nu("2017/01-2017-12"))
    break
