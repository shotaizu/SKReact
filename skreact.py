#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

""" 
    A GUI program to generate information regarding
    reactor neutrinos in Super-Kamiokande: production,
    oscillation, interaction and detection.
"""

from reactor import Reactor
from tkinter import *
from params import *
# Lots of matplotlib gubbins to embed into tkinter
import matplotlib as plt
plt.use("TkAgg") 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import tkinter.ttk as ttk
import pandas as pd
import numpy as np
import os

# Getting all reactor power information into pd df
def extract_reactor_info(react_file_path):
    react_dat = pd.read_excel(react_file_path, header=0)

    # List of reactors
    reactors = []
    for index, data in react_dat.loc[(react_dat["Country Code"] == "JP") 
            | (react_dat["Country Code"] == "KR")].iterrows():
        reactors.append(Reactor(
            data["Country Code"],
            data["Core Name"],
            data["Latitude"],
            data["Longitude"],
            data["Core Type"],
            data["Use Mox?"],
            data["Thermal Power"],
            data[7:]))

    return reactors


def main():
    # Set up tkinter window
    skreact_win = Tk()
    skreact_win.title("SKReact")
    skreact_win.geometry(str(WIN_X) + "x" + str(WIN_Y))

    # Set up the reactor list and names
    reactors = extract_reactor_info(REACT_FILE_PATH)
    reactor_names = [reactor.name for reactor in reactors]

    # Combobox to select reactor to look at
    reactors_lbl = ttk.Label(skreact_win, text = "Reactors (JP or KR)")
    reactors_lbl.grid(column=1,row=0)

    reactors_combo = ttk.Combobox(skreact_win)
    reactors_combo["values"] = reactor_names
    reactors_combo.grid(column=1, row=1)
    reactors_combo.current(0)
    selected_reactor_name = reactors_combo.get()

    # Boxes to select start/end dates
    start_lbl = ttk.Label(skreact_win, text = "Start date:")
    start_lbl.grid(column=3, row=1)
    start_year_combo = ttk.Combobox(skreact_win, width=5)
    start_year_combo["values"] = list(range(2015,2018))
    start_year_combo.current(0)
    start_year_combo.grid(column=4, row=1)
    start_div_lbl = ttk.Label(skreact_win, text = "/")
    start_div_lbl.grid(column=5, row=1)
    start_month_combo = ttk.Combobox(skreact_win, width=2)
    start_month_combo["values"] = list(range(1,13))
    start_month_combo.current(0)
    start_month_combo.grid(column=6, row=1)
    start_year = start_year_combo.get()
    start_month = start_month_combo.get()

    end_lbl = ttk.Label(skreact_win, text = "End date:")
    end_lbl.grid(column=3, row=2)
    end_year_combo = ttk.Combobox(skreact_win, width=5)
    end_year_combo["values"] = list(range(2015,2018))
    end_year_combo.current(1)
    end_year_combo.grid(column=4, row=2)
    end_div_lbl = ttk.Label(skreact_win, text = "/")
    end_div_lbl.grid(column=5, row=2)
    end_month_combo = ttk.Combobox(skreact_win, width=2)
    end_month_combo["values"] = list(range(1,13))
    end_month_combo.current(0)
    end_month_combo.grid(column=6, row=2)
    end_year = end_year_combo.get()
    end_month = end_month_combo.get()

    # Label showing number of nu for period and reactor
    n_nu_title_lbl = ttk.Label(skreact_win, text = "N nu produced in period:")
    n_nu_title_lbl.grid(column=0,row=2)
    n_nu_lbl = ttk.Label(skreact_win, text = "n_nu")
    n_nu_lbl.grid(column=1, row=2)

    # Updating label with n_nu for selected reactor/period
    def update_n_nu(event):
        global start_year
        start_year = start_year_combo.get()
        global start_month
        start_month = start_month_combo.get()
        global end_year
        end_year = end_year_combo.get()
        global end_month
        end_month = end_month_combo.get()
        # print(str(start_year) + "/" + str(start_month))
        # print(str(end_year) + "/" + str(end_month))
        if(end_year < start_year
                or (end_year == start_year and end_month < start_month)):
                    n_nu_lbl["text"] = "Start period after end period"
        else:
            period = str(start_year)+"/"+str(start_month).zfill(2)+"-"+str(end_year)+"/"+str(end_month).zfill(2)
            selected_reactor_name = reactors_combo.get()
            selected_reactor =  next((
                reactor for reactor in reactors if reactor.name == selected_reactor_name), None)
            n_nu = selected_reactor.n_nu(period = period)
            n_nu_lbl['text'] = ("n_nu = %.2E" % n_nu)

    # Binding changing any info to update the number of nu 
    start_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    start_month_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_month_combo.bind("<<ComboboxSelected>>", update_n_nu)
    reactors_combo.bind("<<ComboboxSelected>>", update_n_nu)

    # Plotting the monthly load factors for the selected reactor
    lf_labelframe = ttk.Labelframe(skreact_win, text = "Reactor Monthly Load Factors")
    lf_labelframe.grid(column=1, row=3)
    lf_fig = Figure(figsize=(5,4), dpi=100)
    lf_ax = lf_fig.add_subplot(111)

    lf_monthlys = [reactor.lf_monthly for reactor in reactors]
    lf_monthlys[0].plot(ax=lf_ax)

    lf_canvas = FigureCanvasTkAgg(lf_fig, master=lf_labelframe)
    lf_canvas.draw()
    lf_canvas.get_tk_widget().grid(column=1, row=3)

    # Run the window
    skreact_win.mainloop()


if __name__ == "__main__":
    main()

