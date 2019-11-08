#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
    A GUI program to generate information regarding
    reactor neutrinos in Super-Kamiokande: production,
    oscillation, interaction and detection.
"""

from params import *
from reactor import Reactor
from smear import Smear
from scipy import stats
from tkinter import *
from tkinter import messagebox  # Wildcard doesn't import for some reason
import tkinter.ttk as ttk

# Lots of matplotlib gubbins to embed into tkinter
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.dates as mdates
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import pandas as pd
import numpy as np
from datetime import datetime as dt
import pickle
import random
import time
import math
import copy
import os

# Surpressing a warning bug in numpy library when comparing
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

# For formatting the lf plot later
years = mdates.YearLocator()
months = mdates.MonthLocator(interval=4)
monthsFmt = mdates.DateFormatter("%M")

# The years covered by the files in react_dir
# updated when extract_reactor_info is called
file_year_start = 0
file_year_end = 0

# Getting all reactor power information into pd df
def extract_reactor_info(react_dir):
    # List of reactors, only select JP and KR reactors for now
    reactors = []

    file_names = []
    # Create ordered list of filenames
    for file in os.listdir(react_dir):
        file_name = os.fsdecode(file)
        if file_name.startswith("DB") and file_name.endswith(".xls"):
            file_names.append(os.fsdecode(file))

    # Best to have it in order
    file_names.sort()

    global file_year_start
    global file_year_end
    file_year_start = int(file_names[0][2:6])
    file_year_end = int(file_names[-1][2:6])

    for file_name in file_names:
        # Pull reactor info from first file
        file_year = file_name[2:6]
        print("Importing " + file_name + "...")
        react_dat = pd.read_excel(react_dir + file_name, header=None)
        # Must be in format Country,Name,Lat,Long,Type,Mox?,Pth,LF-monthly
        # Look through xls file's reactors.
        for index, data in react_dat.iterrows():
            # Set the data first to check it makes sense
            # Just reactor info, will deal with erros in lf later
            try:
                data_country = str(data[0]).strip()
                data_name = str(data[1]).strip()
                data_latitude = float(data[2])
                data_longitude = float(data[3])
                data_core_type = str(data[4]).strip()
                data_mox = bool(data[5])
                data_p_th = float(data[6])
                # data_lf = []
                # for month in range(1, 13):
                #     data_lf.append(float(data[6 + month]))
            except ValueError:
                if(VERBOSE_IMPORT_ERR):
                    print("PROBLEM IN .XLS FILE")
                    print("Check file " + file_name + 
                        ", row " + str(index) + " for odd data")
                    input("Press Enter to continue importing...")
                continue
            # if(data_name != "BROWNS FERRY-1"):
            #     continue
            # Check if the reactor on this row is in reactors[]
            in_reactors = False
            for reactor in reactors:
                if reactor.name == data_name:
                    in_reactors = True
                    # Checking if any of the reactor's info has changed
                    # updating to most recent if so
                    if data_country != reactor.country:
                        print("Reactor country has changed in file " +
                            file_name + "!")
                        print("Reactor: " + reactor.name)
                        print(reactor.country + " -> " + data_country)
                        print("Updating...")
                        reactor.country = data_country
                    if data_latitude != reactor.latitude:
                        print("Reactor latitude has changed in file " +
                            file_name + "!")
                        print("Reactor: " + reactor.name)
                        print(str(reactor.latitude) + " -> " + 
                            str(data_latitude))
                        print("Updating...")
                        reactor.latitude = data_latitude
                    if data_longitude != reactor.longitude:
                        print("Reactor longitude has changed in file " +
                            file_name + "!")
                        print("Reactor: " + reactor.name)
                        print(str(reactor.longitude) + " -> " + 
                            str(data_longitude))
                        print("Updating...")
                        reactor.longitude = data_longitude
                    if data_core_type != reactor.core_type:
                        print("Reactor core type has changed in file " +
                            file_name + "!")
                        print("Reactor: " + reactor.name)
                        print(str(reactor.core_type) + " -> " + 
                            str(data_core_type))
                        print("Updating...")
                        reactor.core_type = data_core_type
                    if data_mox != reactor.mox:
                        print("Reactor mox has changed in file " +
                            file_name + "!")
                        print("Reactor: " + reactor.name)
                        print(str(reactor.mox) + " -> " + 
                            str(data_mox))
                        print("Updating...")
                        reactor.mox = data_mox
                    if data_p_th != reactor.p_th:
                        print("Reactor p_th has changed in file " +
                            file_name + "!")
                        print("Reactor: " + reactor.name)
                        print(str(reactor.p_th) + " -> " + 
                            str(data_p_th))
                        print("Updating...")
                        reactor.p_th = data_p_th
                    # Add headers in form used throughout
                    for month in range(1, 13):
                        lf_header = file_year + "/%02i" % month
                        try:
                            # 6 is to skip the reactor data
                            reactor.add_to_lf(lf_header, float(data[6 + month]))
                        except:
                            if(VERBOSE_IMPORT_ERR):
                                print(
                                    "Load factor data for "
                                    + reactor.name
                                    + " in month %i/%02i" % (int(file_year), 
                                    int(month))
                                    + " not float compatible"
                                )
                                print(
                                    "Load factor entry: %s" % data[6 + month]
                                )
                                # input("Press Enter to continue importing...")
                                print("Adding zeros...")
                            reactor.add_to_lf(lf_header, 0.0)

            # Adds reactor if it's not in reactors
            ang_dist = math.sqrt((data_latitude - SK_LAT) ** 2 
                + (data_longitude - SK_LONG) ** 2)
            if not in_reactors and ang_dist < R_THRESH_DEG:
                if(VERBOSE_IMPORT): print("NEW REACTOR: " + data_name) 
                reactors.append(
                    Reactor(
                        data_country,
                        data_name,
                        data_latitude,
                        data_longitude,
                        data_core_type,
                        data_mox,
                        data_p_th,
                        pd.Series([]),  # Load factor
                        default=True,
                        calc_spec=True,
                    )
                )
                # Add up until current file with 0s
                if file_year_start != int(file_year):
                    if(VERBOSE_IMPORT): 
                        print("Retroactively filling data with zeros...")
                for year in range(file_year_start, int(file_year)):
                    for month in range(1, 13):
                        lf_header = "%i/%02i" % (year, month)
                        reactors[-1].add_to_lf(lf_header, 0.0)

                # Now add in current file data
                for month in range(1, 13):
                    lf_header = file_year + "/%02i" % month
                    try:
                        reactors[-1].add_to_lf(lf_header, float(data[6 + month]))
                    except TypeError:
                        if(VERBOSE_IMPORT_ERR):
                            print(
                                "Load factor data for "
                                + reactors[-1].name
                                + " in month %i/%02i" % (int(file_year), month)
                                + " not float compatible"
                            )
                            print(
                                "Load factor entry: %s" % reactors[-1].lf_monthly[lf_header]
                            )
                        reactors[-1].add_to_lf(lf_header, 0.0)
            if not in_reactors and ang_dist >= R_THRESH_DEG:
                if(VERBOSE_IMPORT):
                    print(data[1].strip() + " out of range, skipping...")

        # Checking if reactor isn't present in this file
        # adds zeros for load factor if so
        for reactor in reactors:
            if(reactor.name != "BROWNS FERRY-1"):
                continue
            in_file = False
            for index, data in react_dat.iterrows():
                try:
                    # Check if data[1] is string like
                    data_name = data[1].strip()
                except:
                    # Already handled above
                    continue
                if reactor.name == data_name:
                    in_file = True 
            if not in_file:
                if(VERBOSE_IMPORT):
                    print("NOT IN FILE: " + reactor.name + ", adding zeros")
                for month in range(1, 13):
                    lf_header = file_year + "/%02i" % month
                    reactor.add_to_lf(file_year + "/%02i" % month, 0.0)

        print("...done!")

    reactors.sort(key=lambda x: x.name)

    return reactors


# Creating the list of reactors and buttons
reactors_checkboxes = []
reactors_checkbox_vars = []
reactors_buttons = []
highlighted_reactors = []
highlighted_reactors_names = []

# Making this global so values can be put in file later
reactor_lf_tot = pd.Series()
total_osc_spec = pd.Series()  # Incoming
total_int_spec = pd.Series()  # Interacted
total_wit_spec = pd.Series()  # WIT smeared
highlighted_spec_df = pd.DataFrame()


def main():
    # Try to import geo_nu info
    geo_imported = False
    try:
        # geo_lumi = pd.read_csv(GEO_FILE, sep=" ")
        geo_imported = True
    except FileNotFoundError:
        print("Geo file " + GEO_FILE + " not found!")
        print("Cannot import geoneutrinos information.")

    # Try to calculate smearing matrix
    smear_imported = False
    try:
        wit_smear = Smear(WIT_SMEAR_FILE)
        smear_imported = True
    except FileNotFoundError:
        print("Smear file " + WIT_SMEAR_FILE + " not found!")
        print("Cannot import smearing information.")

    # Extracts from .xls if forced to, does not pickle in this case
    if FORCE_XLS_IMPORT:
        print("Extracting reactor info from " + REACT_DIR)
        default_reactors = extract_reactor_info(REACT_DIR)
    else:
        # Set up the reactor list and names
        try:
            # Pulls from pickle if it exists
            print("Unpickling reactor data from " + REACT_PICKLE + "...")
            with open(REACT_PICKLE, "rb") as pickle_file:
                default_reactors = pickle.load(pickle_file)
            # Have to manually set the whole period from the file
            global file_year_start
            global file_year_end
            file_year_start = int(default_reactors[0].lf_monthly.index[0][:4])
            file_year_end = int(default_reactors[0].lf_monthly.index[-1][:4])
            print("...done!")
        except FileNotFoundError:
            print("Reactor file " + REACT_PICKLE + " not found!")
            print("Extracting reactor info from " + REACT_DIR)
            default_reactors = extract_reactor_info(REACT_DIR)
            with open(REACT_PICKLE, "wb") as pickle_file:
                pickle.dump(default_reactors, pickle_file)

    # Calculate produced spectra for these bins
    for default_reactor in default_reactors:
        default_reactor.set_all_spec()
    default_reactor_names = [reactor.name for reactor in default_reactors]
    reactors = copy.deepcopy(default_reactors)
    reactor_names = default_reactor_names.copy()
    n_reactors = len(reactors)

    # Get oscillation parameters from default (will vary)
    dm_21 = DM_21
    c_13 = C_13_NH
    s_2_12 = S_2_12
    s_13 = S_13_NH

    # INITIALISING ALL MAIN FRAMES
    # =========================================================================

    # Set up tkinter window
    skreact_win = Tk()
    skreact_win.title("SKReact")
    skreact_win.call("tk", "scaling", 1.0)
    # skreact_win.geometry(str(WIN_X) + "x" + str(WIN_Y))
    # A hack to get the OS's name for the default button colour
    test_button = Button()
    default_button_fgc = test_button.cget("fg")


    skreact_title = ttk.Label(
        skreact_win,
        text=(
            "Welcome to SKReact, a GUI reactor neutrino "
            "simulation for Super-Kamiokande"
        ),
    )
    skreact_title.grid(column=0, row=0, columnspan=4)
    title_divider = ttk.Separator(skreact_win, orient=HORIZONTAL)
    title_divider.grid(column=0, row=1, columnspan=5, sticky="ew")

    # map_labelframe = ttk.Labelframe(skreact_win,
    #         text = "Map of SK and Nearby Reactors UNFINISHED")
    # map_labelframe.grid(column=0, row=2)

    # Load factors/ (P/R^2) / event rate etc.
    lf_labelframe = ttk.Labelframe(skreact_win, text="Reactor Monthly Load Factors")
    lf_labelframe.grid(column=0, row=2, sticky=N)

    # List of reactors to select if they contribute
    # Alongside list of buttons to highlight one specifically
    reactors_labelframe = ttk.Labelframe(skreact_win, text="Reactor Selection")
    reactors_labelframe.grid(column=0, row=3, rowspan=1, sticky=N + S + E + W)

    # Incident spectrum.
    int_spec_labelframe = ttk.Labelframe(skreact_win, text="Interaction Spectrum at SK")
    int_spec_labelframe.grid(
        column=1, row=2, columnspan=2, rowspan=2, sticky=N + S + E + W
    )

    # Options for the inc spec, as well as varying osc. params
    int_spec_options_frame = Frame(int_spec_labelframe)
    int_spec_options_frame.grid(column=0, row=1)
    osc_spec_options_labelframe = ttk.Labelframe(
        int_spec_labelframe, text="Vary Osc. Params"
    )
    osc_spec_options_labelframe.grid(column=0, row=2)
    # Ordered flux/n int list of reactors
    reactor_fluxes_labelframe = ttk.LabelFrame(
        int_spec_labelframe, text="Individual Reactor Contributions"
    )
    reactor_fluxes_labelframe.grid(column=1, row=1, rowspan=2, sticky=N + S + E + W)

    # Produced E_spectra
    prod_spec_labelframe = ttk.Labelframe(skreact_win, text="E Spectrum at Production")
    prod_spec_labelframe.grid(column=3, row=2)
    prod_spec_options_labelframe = ttk.Labelframe(
        prod_spec_labelframe, text="View Fuel Contribution"
    )
    prod_spec_options_labelframe.grid(column=0, row=2)

    # Oscillated spectrum.
    osc_spec_labelframe = ttk.Labelframe(skreact_win, text="Oscillated Flux at SK")
    osc_spec_labelframe.grid(column=3, row=3)
    osc_spec_options_frame = Frame(osc_spec_labelframe)
    osc_spec_options_frame.grid(column=0, row=1)

    # =========================================================================

    # Defining the scrollable canvas
    # factor of 25 gives just enough room
    # 30 for extra reactors
    # TODO: Update length dynamically
    reactors_list_canvas = Canvas(
        reactors_labelframe, scrollregion=(0, 0, 400, n_reactors * 30)
    )
    reactors_list_canvas.pack(fill="both", expand=True)

    reactors_scrollbar = Scrollbar(reactors_list_canvas)
    reactors_scrollbar.pack(side=RIGHT, fill=Y)
    reactors_scrollbar.config(command=reactors_list_canvas.yview)

    reactors_list_canvas.config(yscrollcommand=reactors_scrollbar.set)

    # Binding to only scroll canvas if hovering over it
    def _bind_list_to_mousewheel(event):
        reactors_list_canvas.bind_all("<MouseWheel>", _on_mousewheel)

    def _unbind_list_to_mousewheel(event):
        reactors_list_canvas.unbind_all("<MouseWheel>")

    # Dealing with scrolling the reactor list box
    def _on_mousewheel(event):
        reactors_list_canvas.yview_scroll(-1 * (event.delta), "units")

    # Binding scrolling to scroll the reactor list
    # TODO: make it so it only controls it when hovering over
    reactors_list_canvas.bind("<Enter>", _bind_list_to_mousewheel)
    reactors_list_canvas.bind("<Leave>", _unbind_list_to_mousewheel)

    # Select/deselct all reactors in the list then update
    def select_all_reactors(*args):
        for var in reactors_checkbox_vars:
            var.set(1)
        update_n_nu
        return

    def deselect_all_reactors(*args):
        for var in reactors_checkbox_vars:
            var.set(0)
        update_n_nu
        return

    # Buttons to select all, deselect all, add new reactors
    reactor_list_control_frame = Frame(skreact_win)
    select_all_button = Button(text="Select All", command=select_all_reactors)
    select_all_button.grid(in_=reactor_list_control_frame, column=0, row=0)
    deselect_all_button = Button(text="Deselect All", command=deselect_all_reactors)
    deselect_all_button.grid(in_=reactor_list_control_frame, column=1, row=0)

    # Create new generic reactor, add to reactor list, show info
    def add_reactor(*args):
        new_reactor = Reactor(
            "CUSTOM",
            "Custom Reactor",
            35.36,
            138.7,
            "BWR",
            False,
            1000.0,
            pd.Series(0, index=reactors[0].lf_monthly.index),
            False,
            True
        )
        reactors.append(new_reactor)
        create_reactor_list()
        # Button index will alwas be -1 as it was just added
        highlight_reactor(new_reactor, -1)
        show_info(reactors[-1])

    # Adding custom reactors
    add_reactor_button = Button(text="Add Reactor", command=add_reactor)
    add_reactor_button.grid(in_=reactor_list_control_frame, column=2, row=0)

    # Boxes to select start/end dates
    period_labelframe = ttk.Labelframe(skreact_win, text="Period Selection")
    # period_labelframe.pack(in_=reactors_labelframe,side=BOTTOM)
    period_labelframe.grid(in_=lf_labelframe, column=0, row=1)

    # Want this to be above the period selection so needs to pack after
    reactor_list_control_frame.pack(in_=reactors_labelframe, side=BOTTOM)

    start_lbl = ttk.Label(skreact_win, text="From:")
    start_lbl.grid(in_=period_labelframe, column=0, row=0)
    start_year_combo = ttk.Combobox(skreact_win, width=5)
    start_year_combo["values"] = list(range(file_year_start, file_year_end + 1))
    start_year_combo.set(file_year_end)
    start_year_combo.grid(in_=period_labelframe, column=1, row=0)
    start_div_lbl = ttk.Label(skreact_win, text="/")
    start_div_lbl.grid(in_=period_labelframe, column=2, row=0)
    start_month_combo = ttk.Combobox(skreact_win, width=2)
    start_month_combo["values"] = list(range(1, 13))
    start_month_combo.current(0)
    start_month_combo.grid(in_=period_labelframe, column=3, row=0)
    start_year = start_year_combo.get()
    start_month = start_month_combo.get()

    end_lbl = ttk.Label(skreact_win, text=" To: ")
    end_lbl.grid(in_=period_labelframe, column=4, row=0)
    end_year_combo = ttk.Combobox(skreact_win, width=5)
    end_year_combo["values"] = list(range(file_year_start, file_year_end + 1))
    end_year_combo.set(file_year_end)
    end_year_combo.grid(in_=period_labelframe, column=5, row=0)
    end_div_lbl = ttk.Label(skreact_win, text="/")
    end_div_lbl.grid(in_=period_labelframe, column=6, row=0)
    end_month_combo = ttk.Combobox(skreact_win, width=2)
    end_month_combo["values"] = list(range(1, 13))
    end_month_combo.current(11)
    end_month_combo.grid(in_=period_labelframe, column=7, row=0)
    end_year = end_year_combo.get()
    end_month = end_month_combo.get()

    # Listbox of reactor fluxes and names
    reactor_fluxes_scroll = Scrollbar(reactor_fluxes_labelframe)
    reactor_fluxes_scroll.pack(side=RIGHT)
    reactor_fluxes_list = Listbox(reactor_fluxes_labelframe, selectmode="multiple")
    reactor_fluxes_list.pack(side=LEFT, fill=BOTH, expand=1)
    reactor_fluxes_list.config(yscrollcommand=reactor_fluxes_scroll.set)
    reactor_fluxes_scroll.config(command=reactor_fluxes_list.yview)

    # PLOTS ===================================================================
    plt.rc("xtick", labelsize=8)

    # Setting up plot of monthly load factors/power/powerr^2
    # Called lf cause of legacy
    lf_fig = Figure(figsize=(FIG_X, FIG_Y), dpi=100)
    lf_ax = lf_fig.add_subplot(111)
    # Load factor is a %age which occasionally goes over 100
    # lf_ax.set_ylim(0,110)
    lf_tot_ax = lf_ax.twinx()
    lf_canvas = FigureCanvasTkAgg(lf_fig, master=lf_labelframe)
    lf_canvas.get_tk_widget().grid(column=0, row=0)
    lf_options_frame = Frame(lf_labelframe)
    lf_options_frame.grid(column=0, row=2)
    lf_combo = ttk.Combobox(
        lf_options_frame, values=["P/r^2 to SK (MW/km^2)", "P (MW)", "Load Factor (%)"]
    )
    lf_combo.current(0)
    lf_combo.grid(column=0, row=0)
    # lf_toolbar = NavigationToolbar2Tk(lf_canvas, lf_labelframe)

    # Saving the load factor plot
    def save_lf(*args):
        lf_save_win = Toplevel(skreact_win)
        lf_save_win.title("Save LF/P/Pr^-2 Plot")
        filename_label = Label(lf_save_win, text="Filename:")
        filename_label.grid(column=0, row=0)
        filename = Entry(lf_save_win)
        filename.insert(0, "lf_" + time.strftime("%Y%m%d-%H%M%S"))
        filename.grid(column=1, row=0)
        extension = ttk.Combobox(lf_save_win, values=[".pdf", ".png", ".jpg", ".csv"])
        extension.current(0)
        extension.grid(column=2, row=0)

        def save_and_close(*args):
            if extension.get() == ".csv":
                # TODO: Tidy up when OO is implemented
                reactor_lf_tot.to_csv(filename.get() + extension.get())
            else:
                lf_fig.savefig(filename.get() + extension.get())
            # lf_fig.savefig(filename.get() + extension.get())
            lf_save_win.destroy()

        save_button = Button(lf_save_win, text="Save", command=save_and_close)
        save_button.grid(column=0, row=1, columnspan=3)

    # Options to do with the load factor
    # Stack option put in further down after update_n_nu definition
    lf_save_button = Button(lf_options_frame, text="Save as", command=save_lf)
    lf_save_button.grid(column=2, row=0)

    prod_spec_fig = Figure(figsize=(FIG_X, FIG_Y), dpi=100)
    prod_spec_ax = prod_spec_fig.add_subplot(111)
    # osc_spec_ax.set_xlabel("E_nu (MeV)")
    # osc_spec_ax.set_ylabel("n_int (keV^-1 ????)")
    prod_spec_canvas = FigureCanvasTkAgg(prod_spec_fig, master=prod_spec_labelframe)
    # prod_spec_canvas.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1)
    prod_spec_canvas.get_tk_widget().grid(column=0, row=0)
    # prod_spec_toolbar = NavigationToolbar2Tk(prod_spec_canvas, prod_spec_labelframe)
    prod_spec_options_frame = Frame(prod_spec_labelframe)
    prod_spec_options_frame.grid(column=0, row=1)
    prod_spec_label = Label(prod_spec_options_frame, text="N_prod in period = ")
    prod_spec_label.grid(column=2, row=0)

    osc_spec_fig = Figure(figsize=(FIG_X, FIG_Y), dpi=100)
    osc_spec_ax = osc_spec_fig.add_subplot(111)
    osc_spec_canvas = FigureCanvasTkAgg(osc_spec_fig, master=osc_spec_labelframe)
    osc_spec_canvas.get_tk_widget().grid(column=0, row=0)

    int_spec_fig = Figure(figsize=(2 * FIG_X, 2 * FIG_Y), dpi=100)
    int_spec_ax = int_spec_fig.add_subplot(111)
    smear_spec_ax = int_spec_ax
    effs_ax = int_spec_ax.twinx()
    int_spec_canvas = FigureCanvasTkAgg(int_spec_fig, master=int_spec_labelframe)
    int_spec_canvas.get_tk_widget().grid(column=0, row=0, columnspan=2)
    # osc_spec_toolbar = NavigationToolbar2Tk(osc_spec_canvas, osc_spec_labelframe)

    # Generating a nuance file from the oscillated spectrum
    def nuance_osc_spec(*args):
        osc_spec_nuance_win = Toplevel(skreact_win)
        osc_spec_nuance_win.title("Generate POSITRON nuance file for SK")
        filename_label = Label(osc_spec_nuance_win, text="Filename:")
        filename_label.grid(column=0, row=0)
        filename = Entry(osc_spec_nuance_win)
        filename.insert(0, time.strftime("%Y%m%d-%H%M%S"))
        filename.grid(column=1, row=0)
        extension_label = Label(osc_spec_nuance_win, text=".nuance")
        extension_label.grid(column=2, row=0)
        n_events_label = Label(osc_spec_nuance_win, text="n_events:")
        n_events_label.grid(column=0, row=1)
        n_events_entry = Entry(osc_spec_nuance_win)
        n_events_entry.insert(0, "100000")
        n_events_entry.grid(column=1, row=1)

        def nuance_and_close(*args):
            nuance_out = open(filename.get() + ".nuance", "x")
            # Setting up the prob distribution from the spec using rv_discrete
            # rv_discrete only works with integers, so have to map energies to
            # list of integers
            int_map = range(E_BINS)

            # Converting spectrum to probability distribution
            # spec = total_int_spec.to_list()
            # probs = [x/total_int_spec.sum() for x in spec]
            probs = total_int_spec.divide(total_int_spec.sum()).tolist()

            # Set up the dist and generate list from that
            prob_distribution = stats.rv_discrete(
                name="prob_distribution", values=(int_map, probs)
            )
            rvs = prob_distribution.rvs(size=int(n_events_entry.get()))

            for rv in rvs:
                nuance_out.write("begin \n")
                nuance_out.write("info 2 949000 0.0000E+00\n")
                # nuance_out.write("nuance 3 \n")

                theta = random.uniform(0, 2 * math.pi)

                # -1 to get rid of rounding errors causing events
                # to appear nuance_outside the tank
                x = (SK_R - 1) * math.cos(theta)
                y = (SK_R - 1) * math.sin(theta)
                z = random.uniform(-SK_HH, SK_HH)

                vx = random.uniform(-1, 1)
                vy = random.uniform(-1, 1)
                vz = random.uniform(-1, 1)

                px = vx / math.sqrt(vx * vx + vy * vy + vz * vz)
                py = vy / math.sqrt(vx * vx + vy * vy + vz * vz)
                pz = vz / math.sqrt(vx * vx + vy * vy + vz * vz)

                nuance_out.write("vertex %f %f %f 0\n" % (x, y, z))
                nuance_out.write(
                    "track %i %f %f %f %f 0\n"
                    % (POSITRON_PDG, ENERGIES[rv], px, py, pz)
                )
                nuance_out.write("end \n")

            nuance_out.close()
            osc_spec_nuance_win.destroy()

        nuance_button = Button(
            osc_spec_nuance_win, text="Save", command=nuance_and_close
        )
        nuance_button.grid(column=2, row=1)

    # Saving the oscillated spectrum
    def save_osc_spec(*args):
        osc_spec_save_win = Toplevel(skreact_win)
        osc_spec_save_win.title("Save Oscillated Spectrum Plot")
        filename_label = Label(osc_spec_save_win, text="Filename:")
        filename_label.grid(column=0, row=0)
        filename = Entry(osc_spec_save_win)
        filename.insert(0, "osc_" + time.strftime("%Y%m%d-%H%M%S"))
        filename.grid(column=1, row=0)
        extension = ttk.Combobox(
            osc_spec_save_win, values=[".pdf", ".png", ".jpg", ".csv"]
        )
        extension.current(0)
        extension.grid(column=2, row=0)

        def save_and_close(*args):
            if extension.get() == ".csv":
                # TODO: Tidy up when OO is implemented
                total_osc_spec.to_csv(filename.get() + extension.get())
            else:
                osc_spec_fig.savefig(filename.get() + extension.get())
            osc_spec_save_win.destroy()

        save_button = Button(osc_spec_save_win, text="Save", command=save_and_close)
        save_button.grid(column=0, row=1, columnspan=3)

    # Saving the interacted spectrum
    def save_int_spec(*args):
        int_spec_save_win = Toplevel(skreact_win)
        int_spec_save_win.title("Save Incident Spectrum Plot")
        filename_label = Label(int_spec_save_win, text="Filename:")
        filename_label.grid(column=0, row=0)
        filename = Entry(int_spec_save_win)
        filename.insert(0, "int_" + time.strftime("%Y%m%d-%H%M%S"))
        filename.grid(column=1, row=0)
        extension = ttk.Combobox(
            int_spec_save_win, values=[".pdf", ".png", ".jpg", ".csv"]
        )
        extension.current(0)
        extension.grid(column=2, row=0)

        def save_and_close(*args):
            if extension.get() == ".csv":
                # TODO: Tidy up when OO is implemented
                total_int_spec.to_csv(filename.get() + extension.get())
            else:
                int_spec_fig.savefig(filename.get() + extension.get())
            int_spec_save_win.destroy()

        save_button = Button(int_spec_save_win, text="Save", command=save_and_close)
        save_button.grid(column=0, row=1, columnspan=3)

    osc_spec_save_button = Button(
        osc_spec_options_frame, text="Save as", command=save_osc_spec
    )
    osc_spec_save_button.grid(column=2, row=0)
    osc_spec_flx_label = Label(osc_spec_options_frame, text="Total flux in period = ")
    osc_spec_flx_label.grid(column=3, row=0)
    osc_spec_flx_day_label = Label(
        osc_spec_options_frame, text="Avg flux/day in period = "
    )
    osc_spec_flx_day_label.grid(column=3, row=1)
    osc_spec_flx_s_label = Label(osc_spec_options_frame, text="Avg flux/s in period = ")
    osc_spec_flx_s_label.grid(column=3, row=2)

    # Stack option put in further down after update_n_nu definition
    int_spec_save_button = Button(
        int_spec_options_frame, text="Save as", command=save_int_spec
    )
    int_spec_save_button.grid(column=2, row=1)
    int_spec_nuance_button = Button(
        int_spec_options_frame, text="Nuance", command=nuance_osc_spec
    )
    int_spec_nuance_button.grid(column=2, row=2)
    int_spec_int_label = Label(int_spec_options_frame, text="N_int in period = ")
    int_spec_int_label.grid(column=3, row=1)
    int_spec_det_label = Label(int_spec_options_frame, text="N_detected in period = ")
    int_spec_det_label.grid(column=3, row=2)
    # Showing geo_neutrinos NOT YET IMPLEMENTED
    # geo_spec_show_var = IntVar(value=1)
    # geo_spec_show_check = Checkbutton(osc_spec_options_frame,
    #         text="Show Geo",
    #         variable=geo_spec_show_var,
    #         command=update_n_nu)
    # geo_spec_show_check.grid(column=0, row=0)

    # THE MAIN UPDATING FUNCTION
    # =========================================================================
    # =========================================================================
    # TODO: Replot only the lines, not full axes whenever updating
    def update_n_nu(*args):
        update_start = time.time()
        # So it doesn't update before reactor list is set
        try:
            reactors_checkbox_vars[0].get()
        except IndexError:
            return
        start_year = int(start_year_combo.get())
        start_month = int(start_month_combo.get())
        end_year = int(end_year_combo.get())
        end_month = int(end_month_combo.get())
        if end_year < start_year or (
            end_year == start_year and end_month < start_month
        ):
            print("Start period after end period")
        else:
            period = "%i/%02i-%i/%02i" % (start_year, start_month, end_year, end_month)
            # Making datetime objects for calculating average fluxes
            period_start_dt = dt(start_year, start_month, 1)
            # Must be inclusive of last month, so iterate up a month, use day=1
            period_end_dt = dt(
                end_year + ((end_month + 1) // 13), (end_month % 12) + 1, 1
            )
            period_diff_dt = period_end_dt - period_start_dt
            # Clearing old plots an setting labels
            osc_spec_ax.clear()
            int_spec_ax.clear()
            smear_spec_ax.clear()
            effs_ax.clear()
            int_spec_ax.set_xlabel("E_nu [MeV]")
            int_spec_ax.set_ylabel("dN/dE [%g MeV^-1]" % E_INTERVAL)
            osc_spec_ax.set_xlabel("E_nu [MeV]")
            osc_spec_ax.set_ylabel("dN/dE [%g MeV^-1]" % E_INTERVAL)
            prod_spec_ax.clear()
            prod_spec_ax.set_xlabel("E_nu [MeV]")
            prod_spec_ax.set_ylabel("n_prod [MeV^-1 s^-1]")
            lf_ax.clear()
            lf_ax.set_ylabel(lf_combo.get())
            lf_tot_ax.clear()
            # For some reason when plotting dates it uses months as ints
            start_int = (int(start_year) - file_year_start) * 12 + int(start_month) - 1
            end_int = (int(end_year) - file_year_start) * 12 + int(end_month) - 1
            width_int = end_int - start_int

            # LOAD FACTOR PLOTTING
            # =================================================================
            # Make empty series of load factors, sum up for all highlighted
            # Plot load factors from .xls file, may have errors in file to catch
            # Total load factor on same x axis
            # Totals called lf for legacy TODO: change to something more general
            global reactor_lf_tot

            lf_start = time.time()
            reactor_lf_tot = pd.Series(0, index=reactors[0].lf_monthly.index)
            distances = []
            for reactor in reactors:
                distances.append(reactor.dist_to_sk)
                try:
                    # Being explicit with checking combobox values
                    # in case I change them later and don't update this end
                    if lf_combo.get() == "P/r^2 to SK (MW/km^2)":
                        reactor_lf_tot = reactor_lf_tot.add(reactor.p_r_monthly)
                    elif lf_combo.get() == "P (MW)":
                        # if(reactor.p_monthly["2018/04"]<100):
                        #     print(reactor.name)
                        #     print(reactor.p_monthly["2018/04"])
                        reactor_lf_tot = reactor_lf_tot.add(reactor.p_monthly)
                    elif lf_combo.get() == "Load Factor (%)":
                        reactor_lf_tot = reactor_lf_tot.add(reactor.lf_monthly)
                    else:
                        print("ERROR: power/load factor selection not valid")
                        print("Check combobox values in code")
                except TypeError:
                    # Skip over the ones with bad values in the .xls
                    continue

            reactor_lf_tot.plot(ax=lf_tot_ax)
            # print(reactor_lf_tot.loc["2018/02":"2018/06"])

            # To keep the colour same as on osc spec plot where tot is on same ax
            lf_ax.plot(0, 0, alpha=0)

            highlighted_lf_tot = pd.Series(0, index=reactors[0].lf_monthly.index)

            for highlighted_reactor in highlighted_reactors:
                try:
                    if lf_combo.get() == "P/r^2 to SK (MW/km^2)":
                        highlighted_lf = highlighted_reactor.p_r_monthly
                    elif lf_combo.get() == "P (MW)":
                        highlighted_lf = highlighted_reactor.p_monthly
                    elif lf_combo.get() == "Load Factor (%)":
                        highlighted_lf = highlighted_reactor.lf_monthly

                    if lf_stack_var.get():
                        highlighted_lf_tot = highlighted_lf_tot.add(highlighted_lf)
                        highlighted_lf_tot.plot(ax=lf_ax)
                    else:
                        highlighted_lf.plot(ax=lf_ax)
                except TypeError:
                    messagebox.showinfo(
                        "LF Plot Error",
                        "No numeric load factor data to plot! (Check .xls file)"
                        "Reactor: " + highlighted_reactor.name,
                    )

            lf_end = time.time()
            # print("LF runtime = %f" % (lf_end - lf_start))
            # print()

            prod_start = time.time()
            # PRODUCED SPECTRUM PLOTTING
            # =================================================================
            # e_spec on production
            highlighted_e_specs = [
                reactor.prod_spec for reactor in highlighted_reactors
            ]
            spec_errs = [reactor._prod_spec_err() for reactor in highlighted_reactors]
            # Integration
            e_spec_int = 0.0
            # Plotting highlighted fuels
            for highlighted_e_spec, spec_err in zip(highlighted_e_specs, spec_errs):
                for i, fuel in enumerate(highlighted_e_spec.columns.values):
                    # Bit janky relying on order, but same source so fine
                    if plot_fuels_vars[i].get():
                        energies = highlighted_e_spec[fuel].index
                        prod_spec_ax.fill_between(
                            energies,
                            # From when spec_err was errors, not totals
                            # highlighted_e_spec[fuel].add(spec_err[0][fuel]),
                            # highlighted_e_spec[fuel].subtract(spec_err[1][fuel]),
                            spec_err[0][fuel],
                            spec_err[1][fuel],
                            alpha=0.2,
                            color="C%i" % i,
                        )
                        highlighted_e_spec[fuel].plot(ax=prod_spec_ax, color="C%i" % i)
                # Integrating using trap rule
                e_spec_int += np.trapz(
                    highlighted_e_spec["Total"].tolist(), dx=E_INTERVAL
                )

            prod_spec_label["text"] = "N_prod/s @ P_th = %5e" % e_spec_int

            prod_end = time.time()
            # print("Produced runtime = %f" % (prod_end - prod_start))
            # print()

            # INCIDENT SPECTRUM PLOTTING
            # =================================================================
            # Start with empty and add each spectrum
            # This could be done more efficiently
            # TODO: Tidy up once OO
            global total_osc_spec
            global total_int_spec
            global total_wit_spec
            global highlighted_spec_df

            # Add all highlighted spectra to list, concatanate later
            highlighted_osc_specs = []
            highlighted_int_specs = []
            highlighted_colours = []

            # Sum up all spectra
            total_osc_spec = pd.Series(0, index=reactors[0].prod_spec.index)
            if int_spec_offset_var.get() == "e+":
                total_int_spec = pd.Series(0, index=DOWN_ENERGIES)
            else:
                total_int_spec = pd.Series(0, index=ENERGIES)

            # Integration
            spec_start = time.time()
            # print("Spec start...")
            highlight_i = 0

            # Individual reactor total fluxes and names
            reactor_fluxes = []

            for i, reactor in enumerate(reactors):
                if reactors_checkbox_vars[i].get():
                    start = time.time()
                    osc_spec = reactor.osc_spec(
                        dm_21=dm_21_val.get(), s_2_12=s_2_12_val.get(), period=period
                    )
                    end = time.time()
                    # print("Osc spec runtime = %f" % (end-start))

                    start = time.time()
                    int_spec = reactor.int_spec(osc_spec, int_spec_offset_var.get())
                    end = time.time()
                    # print("Int spec runtime = %f" % (end-start))

                    start = time.time()
                    total_osc_spec = total_osc_spec.add(osc_spec)
                    total_int_spec = total_int_spec.add(int_spec)
                    end = time.time()
                    # print("Adding runtime = %f" % (end-start))

                    if reactor in highlighted_reactors:
                        highlighted_osc_specs.append(osc_spec)
                        highlighted_int_specs.append(int_spec)
                        highlighted_colours.append("C%i" % (highlight_i + 1))
                        highlight_i += 1

                    # Add tuple of name and total flux to list
                    reactor_fluxes.append(
                        (
                            reactor.name,
                            np.trapz(osc_spec, dx=E_INTERVAL),
                            np.trapz(int_spec, dx=E_INTERVAL),
                        )
                    )

            # Sort in hi-lo order of fluxes
            reactor_fluxes.sort(key=lambda row: row[1], reverse=True)
            # And put into the listbox (after clearing)
            reactor_fluxes_list.delete(0, END)
            for flux in reactor_fluxes:
                reactor_fluxes_list.insert(
                    END, "%#.4g [cm^-2] | " % (flux[1]) + flux[0]
                )

            spec_end = time.time()
            # print("Total spec runtime = %f" % (spec_end-spec_start))
            # print()

            # Integrating using trap rule
            int_spec_int = np.trapz(total_int_spec.tolist(), dx=E_INTERVAL)
            osc_spec_int = np.trapz(total_osc_spec.tolist(), dx=E_INTERVAL)

            tot_spec_plot_start = time.time()

            # Using C0 so it matches the load factor
            total_osc_spec.plot.area(ax=osc_spec_ax, color="C0", label="Total")
            total_int_spec.plot.area(ax=int_spec_ax, color="C0", label="Total")
            tot_spec_plot_end = time.time()
            # print("Tot plot runtime = %f" % (tot_spec_plot_end-tot_spec_plot_start))
            # print()

            # if(int_spec_eff_var.get()):
            # if(True):
            #     wit_smear.get_effs().plot(
            #         ax = effs_ax,
            #         color = "C3",
            #         label = "Efficiency"
            #     )

            # Exception when nothing is highlighted
            concat_start = time.time()
            try:
                highlighted_osc_spec_df = pd.concat(highlighted_osc_specs, axis=1)
                highlighted_int_spec_df = pd.concat(highlighted_int_specs, axis=1)
                if osc_spec_stack_var.get():
                    highlighted_osc_spec_df.plot.area(
                        ax=osc_spec_ax, color=highlighted_colours
                    )
                else:
                    highlighted_osc_spec_df.plot(
                        ax=osc_spec_ax, color=highlighted_colours
                    )

                if int_spec_stack_var.get():
                    highlighted_int_spec_df.plot.area(
                        ax=int_spec_ax, color=highlighted_colours
                    )
                else:
                    highlighted_int_spec_df.plot(
                        ax=int_spec_ax, color=highlighted_colours
                    )
            except ValueError:
                # Just don't bother concatenating or plotting
                pass
            concat_end = time.time()
            # print("Concat runtime = %f" % (concat_end-concat_start))

            # Plotting smeared spec
            det_spec_int = 0
            if smear_imported:
                smear_spec = wit_smear.smear(total_int_spec, int_spec_offset_var.get())
                smear_spec.plot(ax=smear_spec_ax, color="C3", label="Detected")
                det_spec_int = np.trapz(smear_spec.tolist(), dx=SMEAR_INTERVAL)

            int_spec_int_label["text"] = "N_int in period = %5e" % int_spec_int
            int_spec_det_label["text"] = "N_det in period = %5e" % det_spec_int
            osc_spec_flx_label["text"] = "Total flux in period = %5e" % (osc_spec_int)
            osc_spec_flx_day_label["text"] = "Avg flux/day in period = %5e" % (
                osc_spec_int / period_diff_dt.days
            )
            # For some reason .seconds gives 0 for datetime delta object
            osc_spec_flx_s_label["text"] = "Avg flux/s in period = %5e" % (
                osc_spec_int / period_diff_dt.total_seconds()
            )

            draw_start = time.time()

            # CLEANUP AND DRAWING
            # =================================================================
            prod_spec_ax.legend(loc="lower left")
            prod_spec_ax.set_yscale("log")
            prod_spec_fig.tight_layout()
            prod_spec_canvas.draw()
            # prod_spec_toolbar.update()
            lf_ax.set_ylim(bottom=0)
            lf_ax.set_xlim(start_int - 0.5, end_int + 0.5)
            # lf_ax.xaxis.set_minor_locator(months)
            # lf_ax.xaxis.set_minor_formatter(monthsFmt)
            # lf_fig.autofmt_xdate()
            lf_fig.tight_layout()
            lf_canvas.draw()
            # lf_toolbar.update()
            osc_spec_ax.set_xlim(IBD_MIN, E_MAX)
            osc_spec_ax.set_ylim(bottom=0)
            osc_spec_ax.legend()
            osc_spec_fig.tight_layout()

            int_spec_ax.set_xlim(E_MIN, E_MAX)
            int_spec_ax.set_ylim(bottom=0)
            int_spec_ax.legend(loc="upper right")
            smear_spec_ax.legend(loc="center right")
            int_spec_fig.tight_layout()

            osc_spec_canvas.draw()
            int_spec_canvas.draw()
            # osc_spec_toolbar.update()

            draw_end = time.time()
            # print("Draw runtime = %f" % (draw_end - draw_start))
            # print()

        update_end = time.time()
        # print("Update runtime = %f" % (update_end - update_start))
        # print()
        # print("===========")
        # print()

    # =========================================================================
    # =========================================================================

    # IN PROGRESS
    # Only update load factor plot
    def update_lf(*args):
        pass

    # Only update produced spec plot
    def update_prod(*args):
        pass

    # Only update interacted spec plot
    def update_int(*args):
        pass

    # Update all
    def update_all(*args):
        pass

    # Update load factor and interaction plots
    def update_lf_int(*args):
        pass

    # Choosing whether to stack the oscillation spectra
    int_spec_stack_var = IntVar(value=1)
    int_spec_stack_check = Checkbutton(
        int_spec_options_frame,
        text="Stack",
        variable=int_spec_stack_var,
        command=update_n_nu,
    )
    int_spec_stack_check.grid(column=1, row=0)

    osc_spec_stack_var = IntVar(value=1)
    osc_spec_stack_check = Checkbutton(
        osc_spec_options_frame,
        text="Stack",
        variable=osc_spec_stack_var,
        command=update_n_nu,
    )
    osc_spec_stack_check.grid(column=1, row=0)
    # and if to offset them for nu or e+
    int_spec_offset_var = StringVar()
    int_spec_offset_var.set("e+")
    int_spec_pos_radio = Radiobutton(
        int_spec_options_frame,
        text="e+",
        variable=int_spec_offset_var,
        value="e+",
        command=update_n_nu,
    )
    int_spec_pos_radio.grid(column=2, row=0)
    int_spec_nu_radio = Radiobutton(
        int_spec_options_frame,
        text="nu",
        variable=int_spec_offset_var,
        value="nu",
        command=update_n_nu,
    )
    int_spec_nu_radio.grid(column=3, row=0)

    # Choosing whether to show the efficiency curve
    # int_spec_eff_var = IntVar(value=1)
    # int_spec_eff_check = Checkbutton(int_spec_options_frame,
    #         text="Show effs",
    #         variable=int_spec_eff_var,
    #         command=update_n_nu)
    # int_spec_eff_check.grid(column=1, row=1)

    # Showing (editable) info about a given reactor in new window
    def show_info(reactor):
        reactor_info_win = Toplevel(skreact_win)
        reactor_info_win.title(reactor.name + " Information")
        Label(reactor_info_win, text="Name:").grid(column=0, row=0, sticky=E)
        # For changing name, won't be placed if a def reactor
        name_entry = Entry(reactor_info_win)
        name_entry.insert(0, reactor.name)
        Label(reactor_info_win, text="Country:").grid(column=0, row=1, sticky=E)
        Label(reactor_info_win, text=reactor.country).grid(column=1, row=1, sticky=W)
        # Don't want to change name if default
        # Also can't reset to default if there are no default values
        if reactor.default:
            Label(reactor_info_win, text=reactor.name).grid(column=1, row=0, sticky=W)
        else:
            name_entry.grid(column=1, row=0, sticky=W)

        Label(reactor_info_win, text="Longitude:").grid(column=0, row=2, sticky=E)
        long_entry = Entry(reactor_info_win)
        long_entry.insert(0, reactor.longitude)
        long_entry.grid(column=1, row=2, sticky=W)
        Label(reactor_info_win, text="Latitude:").grid(column=0, row=3, sticky=E)
        lat_entry = Entry(reactor_info_win)
        lat_entry.insert(0, reactor.latitude)
        lat_entry.grid(column=1, row=3, sticky=W)
        Label(reactor_info_win, text="Distance to SK (km)").grid(
            column=0, row=4, sticky=E
        )
        sk_r_entry = Entry(reactor_info_win)
        sk_r_entry.insert(0, "%0.2f" % reactor.dist_to_sk)
        sk_r_entry.grid(column=1, row=4, sticky=W)
        Label(reactor_info_win, text="Core Type?:").grid(column=0, row=5, sticky=E)
        core_type_entry = Entry(reactor_info_win)
        core_type_entry.insert(0, reactor.core_type)
        core_type_entry.grid(column=1, row=5, sticky=W)
        Label(reactor_info_win, text="Uses MOX?:").grid(column=0, row=6, sticky=E)
        mox_check_var = IntVar(value=reactor.mox)
        mox_check = Checkbutton(reactor_info_win, variable=mox_check_var)
        mox_check.grid(column=1, row=6, sticky=W)
        Label(reactor_info_win, text="Thermal Power (Ref/MW):").grid(
            column=0, row=7, sticky=E
        )
        p_th_entry = Entry(reactor_info_win)
        p_th_entry.insert(0, reactor.p_th)
        p_th_entry.grid(column=1, row=7, sticky=W)

        Label(reactor_info_win, text="Monthly Load Factors").grid(column=0, row=8)
        lf_listbox = Listbox(reactor_info_win)

        # Listbox doesn't support row headers/index, so access pd series
        # And combine date and lf into a single string
        for date, lf in reactor.lf_monthly.items():
            lf_listbox.insert(END, date + " - %06.2f" % lf)

        lf_listbox.grid(column=0, row=9)
        lf_entry = Entry(reactor_info_win)
        lf_entry.grid(column=1, row=9)

        # When selecting a listbox item, update the Entry to its value
        def listbox_to_entry(event):
            lf_entry.delete(0, END)
            lf_entry.insert(0, lf_listbox.get(ACTIVE)[-6:])

        # On pressing enter, update listbox selection with Entry value
        def entry_to_listbox(event):
            try:
                item_index = lf_listbox.curselection()[0]
                item_date = lf_listbox.get(ACTIVE)[:7]
                lf_listbox.delete(ACTIVE)
                lf_listbox.insert(
                    item_index, item_date + " - %06.2f" % float(lf_entry.get())
                )
            except IndexError:
                messagebox.showinfo(
                    "LF Input Error", "Please select month to change from list."
                )

        lf_entry.bind("<Return>", entry_to_listbox)

        def lf_series_from_listbox(*args):
            lf_dat = [float(entry[-6:]) for entry in lf_listbox.get(0, END)]
            lf_series = pd.Series(lf_dat, index=reactor.lf_monthly.index)
            return lf_series

        # Updated reactor info with info in the boxes
        def set_reactor_info(*args):
            reactor.set_name(name_entry.get())
            reactor.set_latitude(float(lat_entry.get()))
            reactor.set_longitude(float(long_entry.get()))
            reactor.set_core_type(core_type_entry.get())
            reactor.set_mox(mox_check_var.get())
            reactor.set_p_th(float(p_th_entry.get()))
            reactor.set_lf_monthly(lf_series_from_listbox())
            reactor.set_all_spec()
            # This is definitely needed for custom reactors
            # and I THINK needed otherwise to re-reference in the list
            create_reactor_list()
            update_n_nu()
            return

        # Updates all boxes then sets the reactor info
        def set_reactor_info_def(*args):
            default_reactor = next(
                (x for x in default_reactors if (x.name == reactor.name)), None
            )
            name_entry.delete(0, END)
            name_entry.insert(0, default_reactor.name)
            lat_entry.delete(0, END)
            lat_entry.insert(0, default_reactor.latitude)
            long_entry.delete(0, END)
            long_entry.insert(0, default_reactor.longitude)
            core_type_entry.delete(0, END)
            core_type_entry.insert(0, default_reactor.core_type)
            mox_check_var.set(default_reactor.mox)
            p_th_entry.delete(0, END)
            p_th_entry.insert(0, default_reactor.p_th)
            lf_listbox.delete(0, END)
            for date, lf in default_reactor.lf_monthly.items():
                lf_listbox.insert(END, date + " - %06.2f" % lf)
            set_reactor_info()
            return

        # BROKEN DO NOT USE
        def delete_reactor(*args):
            # Getting index reactor in list
            # This gets completely wrong reactor for some reason
            reactor_i = next(
                (i for i, x in enumerate(reactors) if (x is reactor)), None
            )
            print(i)
            print(reactors[i].name)
            reactors.pop(reactor_i)
            print(reactors[i].name)
            create_reactor_list()
            update_n_nu()
            reactor_info_win.destroy()
            return

        Button(reactor_info_win, text="Update", command=set_reactor_info).grid(
            column=0, row=11
        )
        # Cleaner to re-check down here
        if reactor.default:
            Button(
                reactor_info_win, text="Reset to Def", command=set_reactor_info_def
            ).grid(column=1, row=11)
        # else:
        #     # Give the option to delete custom reactor
        #     Button(reactor_info_win,
        #             text="Delete",
        #             command=delete_reactor
        #             ).grid(column=1,row=11)

    # Creating the list of reactors, once the least of reactors is updated
    def create_reactor_list(*args):

        # First draw map of SK and nearby reactors
        # reac_map_im = plt.imread("japan_map_30_126-43_142.png")
        # map_fig = Figure(figsize=(4,4), dpi=100)
        # map_ax = map_fig.add_subplot(111,label="1")
        # map_ax.axis("off")
        # map_scatter_ax = map_fig.add_subplot(111,label="2")
        # map_canvas = FigureCanvasTkAgg(map_fig,
        #         master=map_labelframe)
        # map_canvas.get_tk_widget().grid(column=0, row=1)
        # map_ax.imshow(reac_map_im)
        reac_lats = [reactor.latitude for reactor in reactors]
        reac_longs = [reactor.longitude for reactor in reactors]
        # Plotting reactor points on top of the image
        # map_scatter_ax.scatter(reac_longs,reac_lats,s=3,c="r")
        # map_scatter_ax.scatter(137.3104,36.4267,s=10,label="Super-Kamiokande")
        # map_scatter_ax.legend()
        # map_scatter_ax.set_xlim(126,142)
        # map_scatter_ax.set_ylim(30,43)
        # map_scatter_ax.patch.set_alpha(0)
        # map_ax.set_xlabel("Latitude (deg)")
        # map_ax.set_ylabel("Longitude (deg)")
        # map_canvas.draw()

        reactors_list_frame = Frame(reactors_list_canvas)
        reactors_list_canvas.create_window(
            (200, 0), window=reactors_list_frame, anchor="n"
        )

        reactors_checkboxes.clear()
        reactors_checkbox_vars.clear()
        reactors_buttons.clear()

        # Header names
        Label(reactors_list_frame, text="Name (Click to Highlight)").grid(
            column=1, row=0, sticky=W
        )
        # Label(reactors_list_frame,text="P_th/MW").grid(column=2,row=0)
        # Label(reactors_list_frame,text="R to SK/km").grid(column=3,row=0)
        # Making the list of reactors and info
        for i, reactor in enumerate(reactors):
            reactors_checkbox_vars.append(IntVar(value=1))
            # reactors_checkbox_vars[i].trace_add("write", update_n_nu)
            reactors_checkboxes.append(
                Checkbutton(
                    reactors_list_frame,
                    variable=reactors_checkbox_vars[i],
                    command=update_n_nu,
                )
            )
            reactors_checkboxes[i].grid(column=0, row=i + 1, sticky=W)
            # Have to explicitly set the index of name to call
            reactors_buttons.append(
                Button(
                    reactors_list_frame,
                    text=reactor.name,
                    command=lambda c=i: highlight_reactor(reactors[c], c),
                )
            )
            reactors_buttons[i].grid(column=1, row=i + 1, sticky=W)
            Button(
                reactors_list_frame,
                text="Info",
                command=lambda c=i: show_info(reactors[c]),
            ).grid(column=2, row=i + 1)

    # Toggled whether the reactor clicked is highlighted or not
    # Taking button index is a hack because this isn't oop
    def highlight_reactor(selected_reactor, button_i):
        # To edit the global, not just create local
        global highlighted_reactors
        global highlighted_reactors_names
        this_button = reactors_buttons[button_i]
        fg_col = this_button.cget("fg")
        # check if it isn't already highlighted
        if fg_col == default_button_fgc:
            this_button.configure(fg="blue")
            highlighted_reactors.append(selected_reactor)
            highlighted_reactors_names.append(selected_reactor.name)
        else:
            this_button.configure(fg=default_button_fgc)
            new_highlighted_reactors = []
            new_highlighted_reactors_names = []
            for reactor in highlighted_reactors:
                if reactor.name != selected_reactor.name:
                    new_highlighted_reactors.append(reactor)
                    new_highlighted_reactors_names.append(reactor.name)
            highlighted_reactors = new_highlighted_reactors.copy()
            highlighted_reactors_names = new_highlighted_reactors_names.copy()

        update_n_nu()

    # Choosing which fuels to show
    plot_fuels_vars = []
    plot_fuels_checks = []
    fuels = reactors[0].prod_spec.columns.values
    for i, fuel in enumerate(fuels):
        # Plot total contribution as default
        if fuel == "Total":
            plot_fuels_vars.append(IntVar(value=1))
        else:
            plot_fuels_vars.append(IntVar(value=0))
        plot_fuels_checks.append(
            Checkbutton(skreact_win, text=fuel, variable=plot_fuels_vars[i])
        )
        plot_fuels_checks[i].grid(in_=prod_spec_options_labelframe, column=i, row=1)

    # And load factors
    lf_stack_var = IntVar(value=1)
    lf_stack_check = Checkbutton(
        lf_options_frame, text="Stack", variable=lf_stack_var, command=update_n_nu
    )
    lf_stack_check.grid(column=1, row=0)

    # Replot when fuels change
    # TODO: Make this plot the fuel itself, save plotting everything each time
    # update_n_nu can then call this
    for i, var in enumerate(plot_fuels_vars):
        plot_fuels_vars[i].trace_add("write", update_n_nu)

    # Reset Osc Params to default
    def reset_osc():
        s_2_12_slider.set(s_2_12)
        dm_21_slider.set(dm_21)
        update_n_nu()

    # Sliders and input to vary the (relevent) osc. params
    # tkinter val shared across slider and input box
    # Have to set values at top so both exist before update is called
    s_2_12_val = DoubleVar(value=s_2_12)
    dm_21_val = DoubleVar(value=dm_21)

    s_2_12_label = ttk.Label(skreact_win, text="Sin^2(theta_12)")
    s_2_12_label.grid(in_=osc_spec_options_labelframe, column=0, row=0)
    s_2_12_val.trace_add("write", update_n_nu)
    s_2_12_slider = Scale(
        skreact_win,
        from_=0,
        to=1,
        resolution=1e-5,
        variable=s_2_12_val,
        orient=HORIZONTAL,
    )
    s_2_12_slider.grid(in_=osc_spec_options_labelframe, column=1, row=0)
    s_2_12_input = Entry(skreact_win, textvariable=s_2_12_val, width=10)
    s_2_12_input.grid(in_=osc_spec_options_labelframe, column=2, row=0)
    dm_21_label = ttk.Label(skreact_win, text="delta m^2_21")
    dm_21_label.grid(in_=osc_spec_options_labelframe, column=0, row=1)
    dm_21_val.trace_add("write", update_n_nu)
    dm_21_slider = Scale(
        skreact_win,
        from_=0,
        to=1e-4,
        resolution=1e-9,
        variable=dm_21_val,
        orient=HORIZONTAL,
    )
    dm_21_slider.grid(in_=osc_spec_options_labelframe, column=1, row=1)
    dm_21_input = Entry(skreact_win, textvariable=dm_21_val, width=10)
    dm_21_input.grid(in_=osc_spec_options_labelframe, column=2, row=1)
    reset_osc_button = Button(skreact_win, text="Reset to default", command=reset_osc)
    reset_osc_button.grid(in_=osc_spec_options_labelframe, column=0, row=2)

    # Binding changing any info to update the number of nu
    lf_combo.bind("<<ComboboxSelected>>", update_n_nu)
    start_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    start_month_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_month_combo.bind("<<ComboboxSelected>>", update_n_nu)

    create_reactor_list()
    # Just to show something on startup
    highlight_reactor(reactors[0], 0)

    update_n_nu()

    # Run the window
    skreact_win.mainloop()


if __name__ == "__main__":
    main()
