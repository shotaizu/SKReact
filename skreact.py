#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    skreact_win.call("tk", "scaling", 1.0)
    # skreact_win.geometry(str(WIN_X) + "x" + str(WIN_Y))

    skreact_title = ttk.Label(skreact_win,
            text = ("Welcome to SKReact, a GUI reactor neutrino "
                "simulation for Super-Kamiokande by Alex Goldsack."))
    skreact_title.grid(column=0, row=0, columnspan=2)
    title_divider = ttk.Separator(skreact_win, orient=HORIZONTAL)
    title_divider.grid(column=0, row=1, columnspan=3, sticky="ew")

    # Set up the reactor list and names
    default_reactors = extract_reactor_info(REACT_FILE_PATH)
    default_reactor_names = [reactor.name for reactor in default_reactors]
    reactors = default_reactors.copy()
    reactor_names = default_reactor_names.copy()
    n_reactors = len(reactors)

    # Get oscillation parameters (will vary)
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
    # 30 for extra reactors, will make it update dynamically at some point
    reactors_list_canvas = Canvas(reactors_labelframe, 
            scrollregion=(0,0,400,n_reactors*30))
    reactors_list_canvas.pack(fill="both", expand=True)

    reactors_scrollbar = Scrollbar(reactors_list_canvas)
    reactors_scrollbar.pack(side=RIGHT, fill=Y)
    reactors_scrollbar.config(command=reactors_list_canvas.yview)

    reactors_list_canvas.config(yscrollcommand=reactors_scrollbar.set)

    # Dealing with scrolling the reactor list box
    def _on_mousewheel(event):
        reactors_list_canvas.yview_scroll(-1*(event.delta), "units")

    # Binding scrolling to scroll the reactor list
    # TODO: make it so it only controls it when hovering over
    reactors_list_canvas.bind_all("<MouseWheel>", _on_mousewheel)

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
    start_year_combo["values"] = list(range(2015,2018))
    start_year_combo.current(0)
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
    end_year_combo["values"] = list(range(2015,2018))
    end_year_combo.current(1)
    end_year_combo.grid(in_=period_labelframe, column=5, row=0)
    end_div_lbl = ttk.Label(skreact_win, text = "/")
    end_div_lbl.grid(in_=period_labelframe, column=6, row=0)
    end_month_combo = ttk.Combobox(skreact_win, width=2)
    end_month_combo["values"] = list(range(1,13))
    end_month_combo.current(0)
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
    map_canvas.draw()

    # Setting up plot of monthly load factors
    lf_labelframe = ttk.Labelframe(skreact_win, 
            text = "Reactor Monthly Load Factors")
    lf_labelframe.grid(column=0, row=3)
    lf_fig = Figure(figsize=(FIG_X,FIG_Y), dpi=100)
    lf_ax = lf_fig.add_subplot(111)
    # Load factor is a %age which occasionally goes over 100
    lf_ax.set_ylim(0,110)
    lf_canvas = FigureCanvasTkAgg(lf_fig, 
            master=lf_labelframe)
    lf_canvas.get_tk_widget().grid(column=2, row=2)

    # And of produced E_spectra
    e_spec_labelframe = ttk.Labelframe(skreact_win, 
            text = "E Spectrum at Production")
    e_spec_labelframe.grid(column=0, row=4)
    e_spec_fig = Figure(figsize=(FIG_X,FIG_Y), dpi=100)
    e_spec_ax = e_spec_fig.add_subplot(111)
    e_spec_canvas = FigureCanvasTkAgg(e_spec_fig, 
            master=e_spec_labelframe)
    e_spec_canvas.get_tk_widget().grid(column=0, row=4)

    # And of oscillated spectrum.
    osc_spec_labelframe = ttk.Labelframe(skreact_win, 
            text = "Prompt E Spectrum at SK")
    osc_spec_labelframe.grid(column=1, row=4)
    osc_spec_fig = Figure(figsize=(FIG_X,FIG_Y), dpi=100)
    osc_spec_ax = osc_spec_fig.add_subplot(111)
    osc_spec_canvas = FigureCanvasTkAgg(osc_spec_fig, 
            master=osc_spec_labelframe)
    osc_spec_canvas.get_tk_widget().grid(column=0, row=6)

    # Updating label with n_nu for highlighted reactor/period
    # TODO: Replot only the lines, not full axes whenever updating
    def update_n_nu(*args):
        global start_year
        start_year = int(start_year_combo.get())
        global start_month
        start_month = int(start_month_combo.get())
        global end_year
        end_year = int(end_year_combo.get())
        global end_month
        end_month = int(end_month_combo.get())
        if(end_year < start_year
                or (end_year == start_year and end_month < start_month)):
                    n_nu_lbl["text"] = "Start period after end period"
        else:
            period = "%i/%02i-%i/%02i" % (start_year, start_month, end_year, end_month)
            highlighted_reactor =  next((
                reactor for reactor in reactors if reactor.name == highlighted_reactor_name), None)
            n_nu = highlighted_reactor.n_nu(period = period)
            n_nu_lbl['text'] = ("n_nu = %.2E" % n_nu)
            osc_spec_ax.clear()
            e_spec_ax.clear()
            lf_ax.clear()
            # For some reason when plotting dates it uses months as ints
            start_int = (int(start_year)-2015)*12+int(start_month)-1
            end_int = (int(end_year)-2015)*12+int(end_month)-1
            width_int = end_int - start_int
            # Box showing period highlighted
            # width to show inclusivity, starting from start of "bin"
            period_box = patches.Rectangle(
                    (start_int-0.5,0), width=width_int+1, height=110, alpha=0.2)
            lf_ax.add_patch(period_box)

            # Plot load factors from .xls file, may have errors in file
            try:
                highlighted_reactor.lf_monthly.plot(ax=lf_ax, marker=".")
            except TypeError:
                messagebox.showinfo("LF Plot Error", 
                        "No numeric load factor data to plot! (Check .xls file)")

            # Plotting incident spectrum
            # Start with empty and add each spectrum
            total_spec = [0]*E_BINS
            for i,reactor in enumerate(reactors):
                if(reactors_checkbox_vars[i]):
                    # print(reactor.name)
                    reactor_spec = reactor.incident_spec(
                        dm_21 = dm_21_val.get(),
                        s_2_12 = s_2_12_val.get()).tolist()
                    total_spec = [f_1 + f_2 for 
                            f_1,f_2 in zip(total_spec,reactor_spec)]
            osc_spec_ax.plot(np.linspace( E_MIN, E_MAX, E_BINS ),total_spec)
            highlighted_reactor.incident_spec(
                    dm_21 = dm_21_val.get(),
                    s_2_12 = s_2_12_val.get()).plot(ax=osc_spec_ax)

            # e_spec on production
            highlighted_e_spec = highlighted_reactor.e_spectra()
            # Plotting highlighted fuels
            for i,fuel in enumerate(highlighted_e_spec.columns.values):
                # Bit janky relying on order, but same source so fine
                if(plot_fuels_vars[i].get()):
                    highlighted_e_spec[fuel].plot(ax=e_spec_ax, color="C%i"%i)
            e_spec_ax.legend(loc="lower left")
            e_spec_ax.set_yscale("log")
            e_spec_canvas.draw()
            lf_fig.autofmt_xdate()
            lf_ax.set_ylim(bottom=0)
            lf_canvas.draw()
            osc_spec_ax.set_xlim(E_MIN,E_MAX)
            osc_spec_ax.set_ylim(bottom=0)
            osc_spec_canvas.draw()

    # Showing (editable) info about a given reactor in new window
    def show_info(reactor):
        reactor_info_win = Toplevel(skreact_win)
        reactor_info_win.title(reactor.name + " Information")
        Label(reactor_info_win,text="Name:").grid(column=0,row=0,sticky=E)
        Label(reactor_info_win,text=reactor.name).grid(column=1,row=0,sticky=W)
        Label(reactor_info_win,text="Country:").grid(column=0,row=1,sticky=E)
        Label(reactor_info_win,text=reactor.country).grid(column=1,row=1,sticky=W)
        Label(reactor_info_win,text="Longitude:").grid(column=0,row=2,sticky=E)
        Label(reactor_info_win,text=reactor.longitude).grid(column=1,row=2,sticky=W)
        Label(reactor_info_win,text="Latitude:").grid(column=0,row=3,sticky=E)
        Label(reactor_info_win,text=reactor.latitude).grid(column=1,row=3,sticky=W)
        Label(reactor_info_win,text="Uses MOX?:").grid(column=0,row=4,sticky=E)
        if(reactor.mox):
            Label(reactor_info_win,text="Yes").grid(column=1,row=4,sticky=W)
        else:
            Label(reactor_info_win,text="No").grid(column=1,row=4,sticky=W)
        Label(reactor_info_win,text="Thermal Power (Ref/MW):").grid(column=0,row=5,sticky=E)
        Label(reactor_info_win,text=reactor.p_th).grid(column=1,row=5,sticky=W)

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
        for i,reactor_name in enumerate(reactor_names):
            reactors_checkbox_vars.append(IntVar(value=1))
            # reactors_checkbox_vars[i].trace_add("write", update_n_nu)
            reactors_checkboxes.append(Checkbutton(reactors_list_frame,
                variable=reactors_checkbox_vars[i],
                command=update_n_nu))
            reactors_checkboxes[i].grid(column=0,row=i+1,sticky=W)
            # Have to explicitly set the index of name to call
            reactors_buttons.append(Button(reactors_list_frame,
                    text=reactor_name,
                    command=lambda c=i: highlight_reactor(reactor_names[c])
                    ))
            reactors_buttons[i].grid(column=1,row=i+1,sticky=W)
            Button(reactors_list_frame,
                    text="Info",
                    command=lambda c=i: show_info(reactors[c])
                    ).grid(column=2,row=i+1)
            # Label(reactors_list_frame,
            #         text="%i"%reactors[i].p_th).grid(column=2,row=i+1)
            # Label(reactors_list_frame,
            #         text="%0.1f"%reactors[i].dist_to_sk).grid(column=3,row=i+1)

    # Highlighting the reactor clicked
    def highlight_reactor(name):
        global highlighted_reactor_name
        highlighted_reactor_name = name
        update_n_nu()

    reactors_checkboxes = []
    reactors_checkbox_vars = []
    reactors_buttons = []
    global highlighted_reactor_name
    highlighted_reactor_name=reactor_names[0]
    create_reactor_list()


    # Choosing which fuels to show
    e_spec_options_labelframe = ttk.Labelframe(skreact_win, 
            text = "View Fuel Contribution")
    e_spec_options_labelframe.grid(column=0, row=6)
    plot_fuels_vars = []
    plot_fuels_checks = []
    fuels = reactors[0].e_spectra().columns.values
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
        plot_fuels_checks[i].grid(in_=e_spec_options_labelframe,column=i,row=1)

    for i,var in enumerate(plot_fuels_vars):
        plot_fuels_vars[i].trace_add("write", update_n_nu)
        # plot_fuels_checks[i].bind("<<CheckbutonSelected>>", update_n_nu)

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
    start_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    start_month_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_year_combo.bind("<<ComboboxSelected>>", update_n_nu)
    end_month_combo.bind("<<ComboboxSelected>>", update_n_nu)

    update_n_nu()

    # Run the window
    skreact_win.mainloop()


if __name__ == "__main__":
    main()

