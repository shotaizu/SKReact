#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
    A GUI program to generate information regarding
    reactor neutrinos in Super-Kamiokande: production,
    oscillation, interaction and detection.
"""

from params import *
from reactor import Reactor
from tkinter import *
from tkinter import messagebox # Wildcard doesn't import for some reason
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
import datetime as dt
import time
import math
import os

# Surpressing a warning bug in numpy library when comparing
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
        react_dat = pd.read_excel(react_dir+file_name, header=None)
        # Must be in format Country,Name,Lat,Long,Type,Mox?,Pth,LF-monthly
        # Add reactor info if first file 
        # Only select Japanese and Korean reactors
        # Too slow otherwise (and changes basically nothing)
        for index, data in react_dat.loc[
                (react_dat[0] == "JP") 
                | (react_dat[0]== "KR")].iterrows():
        # for index, data in react_dat.iterrows():
            # If the reactor on this row is in reactors[]
            added_yet = False
            for reactor in reactors:
                if(reactor.name == data[1].strip()):
                    added_yet = True
                    # Add headers in form used throughout
                    for month in range(1,13):
                        lf_header = file_year + "/%02i" % month
                        # 6 is to skip the reactor data
                        reactor.add_to_lf(lf_header, data[6+month])
                        # reactor.lf_monthly.set_value((file_year + "/%02i" % month),
                        #         data[6+month])
                        try:
                            test_lf_float = float(reactor.lf_monthly[lf_header])
                        except TypeError:
                            print("Load factor data for "
                                    + reactor.name
                                    + " in month %i/%02i" 
                                    % (int(file_year),month)
                                    + " not float compatible")
                            print("Load factor entry: %s" 
                                    % reactor.lf_monthly[lf_header])
                            exit()

            if(not added_yet):
                print("NEW REACTOR: "
                        + data[1].strip())
                reactors.append(Reactor(
                    data[0].strip(), # Country
                    data[1].strip(), # Name
                    data[2], # Lat
                    data[3], # Long
                    data[4], # Type
                    data[5], # Mox?
                    data[6], # Pth
                    pd.Series([]), # Load factor
                    True))
                # Add up until current file with 0s
                if(file_year_start != int(file_year)):
                    print("Retroactively filling data with zeros...")
                for year in range(file_year_start, int(file_year)):
                    for month in range(1,13):
                        lf_header = "%i/%02i" % (year,month)
                        reactors[-1].add_to_lf(lf_header, 0.0)
                        # reactors[-1].lf_monthly.set_value(lf_header,
                        #         0.0)

                # Now add in current file data
                for month in range(1,13):
                    lf_header = file_year + "/%02i" % month
                    reactors[-1].add_to_lf(lf_header, data[6+month])
                    # reactors[-1].lf_monthly.set_value(lf_header,
                    #         data[6+month])
                    try:
                        test_lf_float = float(reactors[-1].lf_monthly[lf_header])
                    except TypeError:
                        print("Load factor data for "
                                + reactors[-1].name
                                + " in month %i/%02i" 
                                % (int(file_year),month)
                                + " not float compatible")
                        print("Load factor entry: %s" 
                                % reactors[-1].lf_monthly[lf_header])
                        exit()
        
        # Checking if reactor isn't present in this file
        # adds zeros for load factor if so
        for reactor in reactors:
            deleted = True
            for index, data in react_dat.loc[
                    (react_dat[0] == "JP") 
                    | (react_dat[0] == "KR")].iterrows():
            # for index, data in react_dat.iterrows():
                if(reactor.name == data[1].strip()):
                    deleted = False
            if(deleted):
                print("NO LONGER TRACKED: "
                + reactor.name
                + ", adding zeros")
                for month in range(1,13):
                    lf_header = file_year + "/%02i" % month
                    reactor.add_to_lf(file_year + "/%02i" % month, 0.0)
                    # reactor.lf_monthly.set_value(file_year + "/%02i" % month,
                    #         0.0)


        print("...done!")

    reactors.sort(key=lambda x: x.name)

    return reactors

# Creating the list of reactors and buttons
reactors_checkboxes = []
reactors_checkbox_vars = []
reactors_buttons = []
highlighted_reactors=[]
highlighted_reactors_names=[]

def main():
    # Set up tkinter window
    skreact_win = Tk()
    skreact_win.title("SKReact")
    skreact_win.call("tk", "scaling", 1.0)
    # skreact_win.geometry(str(WIN_X) + "x" + str(WIN_Y))

    # Try to import geo_nu info
    geo_imported = False
    try:
        geo_lumi = pd.read_csv(GEO_FILE, sep=" ")
        geo_imported = True
    except FileNotFoundError:
        print("File " + GEO_FILE + " not found!")
        print("Cannot import geoneutrinos information.")

    skreact_title = ttk.Label(skreact_win,
            text = ("Welcome to SKReact, a GUI reactor neutrino "
                "simulation for Super-Kamiokande"))
    skreact_title.grid(column=0, row=0, columnspan=2)
    title_divider = ttk.Separator(skreact_win, orient=HORIZONTAL)
    title_divider.grid(column=0, row=1, columnspan=3, sticky="ew")


    # Set up the reactor list and names
    default_reactors = extract_reactor_info(REACT_DIR)
    default_reactor_names = [reactor.name for reactor in default_reactors]
    reactors = default_reactors.copy()
    reactor_names = default_reactor_names.copy()
    n_reactors = len(reactors)

    # Get oscillation parameters from default (will vary)
    dm_21 = DM_21
    c_13 = C_13_NH
    s_2_12 = S_2_12
    s_13 = S_13_NH

    # List of reactors to select if they contribute
    # Alongside list of buttons to highlight one specifically
    reactors_labelframe = ttk.Labelframe(skreact_win, text = "Reactor Selection")
    reactors_labelframe.grid(column=1,row=2,rowspan=2,sticky=N+S+E+W)

    # Defining the scrollable canvas
    # factor of 25 gives just enough room
    # 30 for extra reactors
    # TODO: Update length dynamically
    reactors_list_canvas = Canvas(reactors_labelframe, 
            scrollregion=(0,0,400,n_reactors*30))
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
        reactors_list_canvas.yview_scroll(-1*(event.delta), "units")

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
    reactor_list_control_frame.pack(in_=reactors_labelframe,side=BOTTOM)
    select_all_button = Button(text = "Select All",
            command = select_all_reactors)
    select_all_button.grid(in_=reactor_list_control_frame,column=0,row=0)
    deselect_all_button = Button(text = "Deselect All",
            command = deselect_all_reactors)
    deselect_all_button.grid(in_=reactor_list_control_frame,column=1,row=0)

    # Boxes to select start/end dates
    period_labelframe = ttk.Labelframe(skreact_win, text = "Period Selection")
    period_labelframe.pack(in_=reactors_labelframe,side=BOTTOM)

    start_lbl = ttk.Label(skreact_win, text = "From:")
    start_lbl.grid(in_=period_labelframe, column=0, row=0)
    start_year_combo = ttk.Combobox(skreact_win, width=5)
    start_year_combo["values"] = list(range(file_year_start,file_year_end+1))
    start_year_combo.set(file_year_end)
    start_year_combo.grid(in_=period_labelframe, column=1, row=0)
    start_div_lbl = ttk.Label(skreact_win, text = "/")
    start_div_lbl.grid(in_=period_labelframe, column=2, row=0)
    start_month_combo = ttk.Combobox(skreact_win, width=2)
    start_month_combo["values"] = list(range(1,13))
    start_month_combo.current(0)
    start_month_combo.grid(in_=period_labelframe, column=3, row=0)
    start_year = start_year_combo.get()
    start_month = start_month_combo.get()

    end_lbl = ttk.Label(skreact_win, text = " To: ")
    end_lbl.grid(in_=period_labelframe, column=4, row=0)
    end_year_combo = ttk.Combobox(skreact_win, width=5)
    end_year_combo["values"] = list(range(file_year_start,file_year_end+1))
    end_year_combo.set(file_year_end)
    end_year_combo.grid(in_=period_labelframe, column=5, row=0)
    end_div_lbl = ttk.Label(skreact_win, text = "/")
    end_div_lbl.grid(in_=period_labelframe, column=6, row=0)
    end_month_combo = ttk.Combobox(skreact_win, width=2)
    end_month_combo["values"] = list(range(1,13))
    end_month_combo.current(11)
    end_month_combo.grid(in_=period_labelframe, column=7, row=0)
    end_year = end_year_combo.get()
    end_month = end_month_combo.get()

    # Label showing number of nu for period and reactor
    n_nu_lbl = ttk.Label(skreact_win, text = "n_nu")
    n_nu_lbl.grid(column=0, row=3)

    # PLOTS ===================================================================
    plt.rc('xtick',labelsize=8)

    # Map of SK and nearby reactors
    reac_map_im = plt.imread("japan_map_30_126-43_142.png")
    map_labelframe = ttk.Labelframe(skreact_win, 
            text = "Map of SK and Nearby Reactors UNFINISHED")
    map_labelframe.grid(column=0, row=2)
    map_fig = Figure(figsize=(4,4), dpi=100)
    map_ax = map_fig.add_subplot(111,label="1")
    map_ax.axis("off")
    map_scatter_ax = map_fig.add_subplot(111,label="2")
    map_canvas = FigureCanvasTkAgg(map_fig, 
            master=map_labelframe)
    map_canvas.get_tk_widget().grid(column=0, row=1)
    map_ax.imshow(reac_map_im)
    reac_lats = [reactor.latitude for reactor in reactors]
    reac_longs = [reactor.longitude for reactor in reactors]
    # Plotting reactor points on top of the image
    map_scatter_ax.scatter(reac_longs,reac_lats,s=3,c="r")
    map_scatter_ax.scatter(137.3104,36.4267,s=10,label="Super-Kamiokande")
    map_scatter_ax.legend()
    map_scatter_ax.set_xlim(126,142)
    map_scatter_ax.set_ylim(30,43)
    map_scatter_ax.patch.set_alpha(0)
    map_ax.set_xlabel("Latitude (deg)")
    map_ax.set_ylabel("Longitude (deg)")
    map_canvas.draw()

    # Setting up plot of monthly load factors/power/powerr^2
    # Called lf cause of legacy
    lf_labelframe = ttk.Labelframe(skreact_win, 
            text = "Reactor Monthly Load Factors")
    lf_labelframe.grid(column=0, row=3)
    lf_fig = Figure(figsize=(FIG_X,FIG_Y), dpi=100)
    lf_ax = lf_fig.add_subplot(111)
    # Load factor is a %age which occasionally goes over 100
    # lf_ax.set_ylim(0,110)
    lf_tot_ax = lf_ax.twinx()
    lf_canvas = FigureCanvasTkAgg(lf_fig, 
            master=lf_labelframe)
    lf_canvas.get_tk_widget().grid(column=0, row=0)
    lf_options_frame = Frame(lf_labelframe)
    lf_options_frame.grid(column=0, row=1)
    lf_combo = ttk.Combobox(
            lf_options_frame,
            values=[
                "P/r^2 to SK (MW/km^2)",
                "P (MW)",
                "Load Factor (%)"])
    lf_combo.current(0)
    lf_combo.grid(column=0,row=0)
    # lf_toolbar = NavigationToolbar2Tk(lf_canvas, lf_labelframe)

    # Saving the load factor plot 
    def save_lf(*args):
        lf_save_win = Toplevel(skreact_win)
        lf_save_win.title("Save LF/P/Pr^-2 Plot")
        filename_label = Label(
                lf_save_win,
                text = "Filename:")
        filename_label.grid(column=0,row=0)
        filename = Entry(lf_save_win)
        filename.insert(0,"lf_" + time.strftime("%Y%m%d-%H%M%S"))
        filename.grid(column=1,row=0)
        extension = ttk.Combobox(
                lf_save_win,
                values = [
                    ".pdf",
                    ".png",
                    ".jpg"])
        extension.current(0)
        extension.grid(column=2,row=0)
        def save_and_close(*args):
            lf_fig.savefig(filename.get() + extension.get())
            lf_save_win.destroy()
                    
        save_button = Button(
                lf_save_win,
                text="Save",
                command = save_and_close
                )
        save_button.grid(column=0,row=1,columnspan=3)

    # Options to do with the load factor 
    # Stack option put in further down after update_n_nu definition
    lf_save_button = Button(lf_options_frame,
            text = "Save as", 
            command=save_lf)
    lf_save_button.grid(column=2,row=0)

    # And of produced E_spectra
    prod_spec_labelframe = ttk.Labelframe(skreact_win, 
            text = "E Spectrum at Production")
    prod_spec_labelframe.grid(column=0, row=4)
    prod_spec_fig = Figure(figsize=(FIG_X,FIG_Y), dpi=100)
    prod_spec_ax = prod_spec_fig.add_subplot(111)
    # osc_spec_ax.set_xlabel("E_nu (MeV)")
    # osc_spec_ax.set_ylabel("n_int (keV^-1 ????)")
    prod_spec_canvas = FigureCanvasTkAgg(prod_spec_fig, 
            master=prod_spec_labelframe)
    # prod_spec_canvas.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1)
    prod_spec_canvas.get_tk_widget().grid(column=0, row=0)
    # prod_spec_toolbar = NavigationToolbar2Tk(prod_spec_canvas, prod_spec_labelframe)
    prod_spec_options_frame = Frame(prod_spec_labelframe)
    prod_spec_options_frame.grid(column=0, row=1)
    prod_spec_label = Label(prod_spec_options_frame,
            text = "N_int in period = ")
    prod_spec_label.grid(column=2,row=0)

    # And of oscillated spectrum.
    osc_spec_labelframe = ttk.Labelframe(skreact_win, 
            text = "Interaction E Spectrum at SK")
    osc_spec_labelframe.grid(column=1, row=4)
    osc_spec_fig = Figure(figsize=(FIG_X,FIG_Y), dpi=100)
    osc_spec_ax = osc_spec_fig.add_subplot(111)
    osc_spec_canvas = FigureCanvasTkAgg(osc_spec_fig, 
            master=osc_spec_labelframe)
    osc_spec_canvas.get_tk_widget().grid(column=0, row=0)
    # osc_spec_toolbar = NavigationToolbar2Tk(osc_spec_canvas, osc_spec_labelframe)

    # Saving the oscillated spectrum
    def save_osc_spec(*args):
        osc_spec_save_win = Toplevel(skreact_win)
        osc_spec_save_win.title("Save Incident Spectrum Plot")
        filename_label = Label(
                osc_spec_save_win,
                text = "Filename:")
        filename_label.grid(column=0,row=0)
        filename = Entry(osc_spec_save_win)
        filename.insert(0,"osc_" + time.strftime("%Y%m%d-%H%M%S"))
        filename.grid(column=1,row=0)
        extension = ttk.Combobox(
                osc_spec_save_win,
                values = [
                    ".pdf",
                    ".png",
                    ".jpg"])
        extension.current(0)
        extension.grid(column=2,row=0)
        def save_and_close(*args):
            osc_spec_fig.savefig(filename.get() + extension.get())
            osc_spec_save_win.destroy()
                    
        save_button = Button(
                osc_spec_save_win,
                text="Save",
                command = save_and_close
                )
        save_button.grid(column=0,row=1,columnspan=3)

    # Options to do with the incident spectrum
    # Stack option put in further down after update_n_nu definition
    osc_spec_options_frame = Frame(osc_spec_labelframe)
    osc_spec_options_frame.grid(column=0, row=1)
    osc_spec_save_button = Button(osc_spec_options_frame,
            text = "Save as", 
            command=save_osc_spec)
    osc_spec_save_button.grid(column=2,row=0)
    osc_spec_int_label = Label(osc_spec_options_frame,
            text = "N_int in period = ")
    osc_spec_int_label.grid(column=3,row=0)

    # THE MAIN UPDATING FUNCTION
    # =========================================================================
    # =========================================================================
    # TODO: Replot only the lines, not full axes whenever updating
    def update_n_nu(*args):
        # So it doesn't update before reactor list is set
        try:
            reactors_checkbox_vars[0].get()
        except IndexError:
            return
        start_year = int(start_year_combo.get())
        start_month = int(start_month_combo.get())
        end_year = int(end_year_combo.get())
        end_month = int(end_month_combo.get())
        if(end_year < start_year
                or (end_year == start_year and end_month < start_month)):
                    n_nu_lbl["text"] = "Start period after end period"
        else:
            period = "%i/%02i-%i/%02i" % (start_year, start_month, end_year, end_month)
            # n_nu = highlighted_reactor.n_nu(period = period)
            # n_nu_lbl['text'] = ("n_nu = %.2E" % n_nu)
            # Clearing old plots an setting labels
            osc_spec_ax.clear()
            osc_spec_ax.set_xlabel("E_nu [MeV]")
            osc_spec_ax.set_ylabel("n_int [MeV^-1]")
            prod_spec_ax.clear()
            prod_spec_ax.set_xlabel("E_nu [MeV]")
            prod_spec_ax.set_ylabel("n_prod [MeV^-1 s^-1]")
            lf_ax.clear()
            lf_ax.set_ylabel(lf_combo.get())
            lf_tot_ax.clear()
            # For some reason when plotting dates it uses months as ints
            start_int = (int(start_year)-file_year_start)*12+int(start_month)-1
            end_int = (int(end_year)-file_year_start)*12+int(end_month)-1
            width_int = end_int - start_int

            # LOAD FACTOR PLOTTING
            # =================================================================
            # Make empty series of load factors, sum up for all highlighted
            # Plot load factors from .xls file, may have errors in file to catch
            # Total load factor on same x axis
            # Totals called lf for legacy TODO: change to something more general
            reactor_lf_tot = pd.Series(0,
                    index=reactors[0].lf_monthly.index)
            for reactor in reactors:
                try:
                    # Being explicit with checking combobox values
                    # in case I change them later and don't update this end
                    if(lf_combo.get() == "P/r^2 to SK (MW/km^2)"):
                        reactor_lf_tot = reactor_lf_tot.add(reactor.p_r_monthly)
                    elif(lf_combo.get() == "P (MW)"):
                        reactor_lf_tot = reactor_lf_tot.add(reactor.p_monthly)
                    elif(lf_combo.get() == "Load Factor (%)"):
                        reactor_lf_tot = reactor_lf_tot.add(reactor.lf_monthly)
                    else:
                        print("ERROR: power/load factor selection not valid")
                        print("Check combobox values in code")
                except TypeError:
                    # Skip over the ones with bad values in the .xls
                    continue

            reactor_lf_tot.plot(ax=lf_tot_ax)

            # Box showing period highlighted
            # width to show inclusivity, starting from start of "bin"
            period_box = patches.Rectangle(
                    (start_int-0.5,0), 
                    width=width_int+1, 
                    height=reactor_lf_tot.max(), 
                    alpha=0.2)
            lf_tot_ax.add_patch(period_box)

            # To keep the colour same as on osc spec plot where tot is on same ax
            lf_ax.plot(0,0,alpha=0) 
            
            highlighted_lf_tot = pd.Series(0, 
                    index=reactors[0].lf_monthly.index)

            for highlighted_reactor in highlighted_reactors:
                try:
                    if(lf_combo.get() == "P/r^2 to SK (MW/km^2)"):
                        highlighted_lf = highlighted_reactor.p_r_monthly
                    elif(lf_combo.get() == "P (MW)"):
                        highlighted_lf = highlighted_reactor.p_monthly
                    elif(lf_combo.get() == "Load Factor (%)"):
                        highlighted_lf = highlighted_reactor.lf_monthly

                    if(lf_stack_var.get()):
                        highlighted_lf_tot = highlighted_lf_tot.add(
                                highlighted_lf)
                        highlighted_lf_tot.plot(ax=lf_ax)
                    else:
                        highlighted_lf.plot(ax=lf_ax)
                except TypeError:
                    messagebox.showinfo("LF Plot Error", 
                            "No numeric load factor data to plot! (Check .xls file)"
                            "Reactor: " + highlighted_reactor.name)

            # PRODUCED SPECTRUM PLOTTING
            # =================================================================
            # e_spec on production
            highlighted_e_specs = [reactor.e_spectra
                    for reactor in highlighted_reactors]
            # Integration
            e_spec_int = 0.0
            # Plotting highlighted fuels
            for highlighted_e_spec in highlighted_e_specs:
                for i,fuel in enumerate(highlighted_e_spec.columns.values):
                    # Bit janky relying on order, but same source so fine
                    if(plot_fuels_vars[i].get()):
                        highlighted_e_spec[fuel].plot(ax=prod_spec_ax, color="C%i"%i)
                # Integrating using trap rule
                e_spec_int += np.trapz(highlighted_e_spec["Total"].tolist(),
                    dx = E_INTERVAL)

            prod_spec_label["text"] = "N_prod/s @ P_th = %5e" % e_spec_int 

            # INCIDENT SPECTRUM PLOTTING
            # =================================================================
            # Start with empty and add each spectrum
            # This could be done more efficiently
            # TODO: Change to pd series
            total_spec = pd.Series(0, 
                    index = reactors[0].e_spectra.index)
            stacked_highlighted_spec = pd.Series(0,
                    index = reactors[0].e_spectra.index)
            # Highlighted plot colour index
            c_i = 1
            # Integration
            int_spec_int = 0.0
            for i,reactor in enumerate(reactors):
                if(reactors_checkbox_vars[i].get()):
                    osc_spec = reactor.oscillated_spec(
                        dm_21 = dm_21_val.get(),
                        s_2_12 = s_2_12_val.get(),
                        period = period)
                    reactor_spec = reactor.incident_spec(osc_spec)

                    total_spec = total_spec.add(reactor_spec)

                    # Integrating using trap rule
                    int_spec_int += np.trapz(reactor_spec.tolist(),
                        dx = E_INTERVAL)
                if(reactor.name in highlighted_reactors_names):
                    highlighted_osc_spec = reactor.oscillated_spec(
                            dm_21 = dm_21_val.get(),
                            s_2_12 = s_2_12_val.get(),
                            period = period)
                    highlighted_spec = reactor.incident_spec(highlighted_osc_spec)
                    if(osc_spec_stack_var.get()):
                        stacked_highlighted_spec = stacked_highlighted_spec.add(
                                highlighted_spec)
                        stacked_highlighted_spec.plot(ax = osc_spec_ax,
                                label = reactor.name,
                                color = "C%i"%c_i)
                    else:
                        highlighted_spec.plot(ax = osc_spec_ax,
                                label = reactor.name,
                                color = "C%i"%i)
                    c_i+=1

            # Using C0 so it matches the load factor
            total_spec.plot(ax = osc_spec_ax, color = "C0", label = "Total")
            
            osc_spec_int_label["text"] = "N_int in period = %5e" % int_spec_int 

            # CLEANUP AND DRAWING 
            # =================================================================
            prod_spec_ax.legend(loc="lower left")
            prod_spec_ax.set_yscale("log")
            # prod_spec_fig.tight_layout()
            prod_spec_canvas.draw()
            # prod_spec_toolbar.update()
            # lf_ax.xaxis.set_major_locator(months)
            # lf_ax.xaxis.set_major_formatter(monthsFmt)
            lf_fig.autofmt_xdate()
            # lf_fig.tight_layout()
            lf_ax.set_ylim(bottom=0)
            lf_canvas.draw()
            # lf_toolbar.update()
            osc_spec_ax.set_xlim(E_MIN,E_MAX)
            osc_spec_ax.set_ylim(bottom=0)
            osc_spec_ax.legend()
            osc_spec_fig.tight_layout()
            osc_spec_canvas.draw()
            # osc_spec_toolbar.update()

    # =========================================================================
    # =========================================================================


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


    # This can go in the show_info function maybe? 
    # def set_reactor_info(name,
    #         latitude,
    #         longitude,
    #         core_type,
    #         mox,
    #         p_th):
    def set_reactor_info(*args):
        return

    # Showing (editable) info about a given reactor in new window
    def show_info(reactor):
        reactor_info_win = Toplevel(skreact_win)
        reactor_info_win.title(reactor.name + " Information")
        Label(reactor_info_win,text="Name:").grid(column=0,row=0,sticky=E)
        Label(reactor_info_win,text=reactor.name).grid(column=1,row=0,sticky=W)
        Label(reactor_info_win,text="Country:").grid(column=0,row=1,sticky=E)
        Label(reactor_info_win,text=reactor.country).grid(column=1,row=1,sticky=W)
        Label(reactor_info_win,text="Longitude:").grid(column=0,row=2,sticky=E)
        long_entry = Entry(reactor_info_win)
        long_entry.insert(0, reactor.longitude)
        long_entry.grid(column=1,row=2,sticky=W)
        Label(reactor_info_win,text="Latitude:").grid(column=0,row=3,sticky=E)
        lat_entry = Entry(reactor_info_win)
        lat_entry.insert(0, reactor.latitude)
        lat_entry.grid(column=1,row=3,sticky=W)
        Label(reactor_info_win,text="Distance to SK (km)").grid(column=0,row=4,sticky=E)
        sk_r_entry = Entry(reactor_info_win)
        sk_r_entry.insert(0, "%0.2f"%reactor.dist_to_sk)
        sk_r_entry.grid(column=1,row=4,sticky=W)
        Label(reactor_info_win,text="Uses MOX?:").grid(column=0,row=5,sticky=E)
        mox_check_var = IntVar(value=reactor.mox)
        mox_check = Checkbutton(reactor_info_win,variable=mox_check_var)
        mox_check.grid(column=1,row=5,sticky=W)
        Label(reactor_info_win,text="Thermal Power (Ref/MW):").grid(column=0,row=6,sticky=E)
        p_th_entry = Entry(reactor_info_win)
        p_th_entry.insert(0, reactor.p_th)
        p_th_entry.grid(column=1,row=6,sticky=W)

        Label(reactor_info_win, text="Monthly Load Factors").grid(column=0,row=7)
        lf_listbox = Listbox(reactor_info_win)

        # Listbox doesn't support row headers/index, so access pd series
        # And combine date and lf into a single string
        for date, lf in reactor.lf_monthly.items():
            lf_listbox.insert(END, date + " - %06.2f"%lf)

        lf_listbox.grid(column=0,row=8)
        lf_entry = Entry(reactor_info_win)
        lf_entry.grid(column=1,row=8)

        # When selecting a listbox item, update the Entry to its value
        def listbox_to_entry(event):
            lf_entry.delete(0, END)
            lf_entry.insert(0, lf_listbox.get(ACTIVE)[-6:])

        # On pressing enter, update listbox selection with Entry value
        def entry_to_listbox(event):
            try:
                item_index = lf_listbox.curselection()[0]
            except IndexError:
                messagebox.showinfo("LF Input Error",
                        "Please select month to change from list.")
            item_date = lf_listbox.get(ACTIVE)[:7]
            lf_listbox.delete(ACTIVE)
            lf_listbox.insert(item_index, 
                    item_date + " - %06.2f"%float(lf_entry.get()))

        lf_entry.bind("<Return>", entry_to_listbox)

        Button(reactor_info_win,
                text="Update (NOT IMPLEMENTED)",
                command=set_reactor_info
                ).grid(column=0,row=10)
        Button(reactor_info_win,
                text="Reset to Def",
                command=set_reactor_info
                ).grid(column=1,row=10)

    # Creating the list of reactors, once the least of reactors is updated
    def create_reactor_list(*args):
        reactors_list_frame = Frame(reactors_list_canvas)
        reactors_list_canvas.create_window((200,0), window=reactors_list_frame, anchor="n")

        reactors_checkboxes.clear()
        reactors_checkbox_vars.clear()
        reactors_buttons.clear()

        # Header names
        Label(reactors_list_frame,text="Name (Click to Highlight)").grid(column=1,row=0,sticky=W)
        # Label(reactors_list_frame,text="P_th/MW").grid(column=2,row=0)
        # Label(reactors_list_frame,text="R to SK/km").grid(column=3,row=0)
        # Making the list of reactors and info
        for i,reactor in enumerate(reactors):
            reactors_checkbox_vars.append(IntVar(value=1))
            # reactors_checkbox_vars[i].trace_add("write", update_n_nu)
            reactors_checkboxes.append(Checkbutton(reactors_list_frame,
                variable=reactors_checkbox_vars[i],
                command=update_n_nu))
            reactors_checkboxes[i].grid(column=0,row=i+1,sticky=W)
            # Have to explicitly set the index of name to call
            reactors_buttons.append(Button(reactors_list_frame,
                    text=reactor.name,
                    command=lambda c=i: highlight_reactor(reactors[c],c)
                    ))
            reactors_buttons[i].grid(column=1,row=i+1,sticky=W)
            Button(reactors_list_frame,
                    text="Info",
                    command=lambda c=i: show_info(reactors[c])
                    ).grid(column=2,row=i+1)

    # Toggled whether the reactor clicked is highlighted or not
    # Taking button index is a hack because this isn't oop
    def highlight_reactor(selected_reactor,button_i):
        # To edit the global, not just create local
        global highlighted_reactors
        global highlighted_reactors_names
        this_button = reactors_buttons[button_i]
        fg_col = this_button.cget("fg")
        # check if it isn't already highlighted
        if(fg_col == "systemButtonText"):
            this_button.configure(fg = "blue") 
            highlighted_reactors.append(selected_reactor)
            highlighted_reactors_names.append(selected_reactor.name)
        else:
            this_button.configure(fg = "systemButtonText")
            new_highlighted_reactors = []
            new_highlighted_reactors_names = []
            for reactor in highlighted_reactors:
                if(reactor.name != selected_reactor.name):
                    new_highlighted_reactors.append(reactor)
                    new_highlighted_reactors_names.append(reactor.name)
            highlighted_reactors = new_highlighted_reactors.copy()
            highlighted_reactors_names = new_highlighted_reactors_names.copy()
        
        update_n_nu()


    # Choosing which fuels to show
    prod_spec_options_labelframe = ttk.Labelframe(skreact_win, 
            text = "View Fuel Contribution")
    prod_spec_options_labelframe.grid(column=0, row=6)
    plot_fuels_vars = []
    plot_fuels_checks = []
    fuels = reactors[0].e_spectra.columns.values
    for i,fuel in enumerate(fuels):
        # Plot total contribution as default
        if(fuel == "Total"):
            plot_fuels_vars.append(IntVar(value=1))
        else:
            plot_fuels_vars.append(IntVar(value=0))
        plot_fuels_checks.append(
                Checkbutton(skreact_win, text=fuel, 
                    variable=plot_fuels_vars[i],
                    ))
        plot_fuels_checks[i].grid(in_=prod_spec_options_labelframe,column=i,row=1)

    # Choosing whether to stack the oscillation spectra
    osc_spec_stack_var = IntVar(value=1)
    osc_spec_stack_check = Checkbutton(osc_spec_options_frame, 
            text="Stack",
            variable=osc_spec_stack_var,
            command=update_n_nu)
    osc_spec_stack_check.grid(column=1, row=0)
    # Showing geo_neutrinos NOT YET IMPLEMENTED
    geo_spec_show_var = IntVar(value=1)
    geo_spec_show_check = Checkbutton(osc_spec_options_frame, 
            text="Show Geo",
            variable=geo_spec_show_var,
            command=update_n_nu)
    geo_spec_show_check.grid(column=0, row=0)

    # And load factors 
    lf_stack_var = IntVar(value=1)
    lf_stack_check = Checkbutton(lf_options_frame, 
            text="Stack",
            variable=lf_stack_var,
            command=update_n_nu)
    lf_stack_check.grid(column=1, row=0)

    # Replot when fuels change
    # TODO: Make this plot the fuel itself, save plotting everything each time
    # update_n_nu can then call this
    for i,var in enumerate(plot_fuels_vars):
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

    osc_spec_options_labelframe = ttk.Labelframe(skreact_win, 
            text = "Vary Osc. Params")
    osc_spec_options_labelframe.grid(column=1,row=6)
    s_2_12_label = ttk.Label(skreact_win, text = "Sin^2(theta_12)")
    s_2_12_label.grid(in_=osc_spec_options_labelframe,column=0,row=0)
    s_2_12_val.trace_add("write", update_n_nu)
    s_2_12_slider = Scale(skreact_win, 
            from_=0, 
            to=1, 
            resolution=1e-5,
            variable=s_2_12_val,
            orient=HORIZONTAL)
    s_2_12_slider.grid(in_=osc_spec_options_labelframe,column=1,row=0)
    s_2_12_input = Entry(skreact_win,textvariable=s_2_12_val,width=10)
    s_2_12_input.grid(in_=osc_spec_options_labelframe,column=2,row=0)
    dm_21_label = ttk.Label(skreact_win, text = "delta m^2_21")
    dm_21_label.grid(in_=osc_spec_options_labelframe,column=0,row=1)
    dm_21_val.trace_add("write", update_n_nu)
    dm_21_slider = Scale(skreact_win, 
            from_=0, 
            to=1e-4, 
            resolution=1e-9,
            variable=dm_21_val,
            orient=HORIZONTAL)
    dm_21_slider.grid(in_=osc_spec_options_labelframe,column=1,row=1)
    dm_21_input = Entry(skreact_win,textvariable=dm_21_val,width=10)
    dm_21_input.grid(in_=osc_spec_options_labelframe,column=2,row=1)
    reset_osc_button = Button(skreact_win,
            text = "Reset to default",
            command = reset_osc)
    reset_osc_button.grid(in_=osc_spec_options_labelframe,column=0,row=2)

    # Binding changing any info to update the number of nu 
    lf_combo.bind("<<ComboboxSelected>>", update_n_nu)
    start_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    start_month_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_month_combo.bind("<<ComboboxSelected>>", update_n_nu)

    create_reactor_list()
    # Just to show something on startup
    highlight_reactor(reactors[0],0)

    update_n_nu()

    # Run the window
    skreact_win.mainloop()


if __name__ == "__main__":
    main()

