from params import *
from reactor import Reactor
from smear import Smear
from tkinter import *
from tkinter import filedialog
import tkinter.ttk as ttk
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

# Window to handle data import and fitting
def fit_win(import_filename,reactors,period,wit_smear):
    try:
        import_dat_df = pd.read_csv(import_filename)
    except:
        print("Cannot read import data!")
        return

    # Change to series
    import_dat_df.columns = ["energy","bin_content"]
    import_dat = import_dat_df.set_index("energy")["bin_content"]

    # Area normalise
    import_dat_int = np.trapz(import_dat)
    import_dat_norm = import_dat.apply(
        lambda x: x/import_dat_int)

    import_dat_norm_int = np.trapz(import_dat_norm)

    # Now make window
    # fit_win = Toplevel(skreact_win)
    fit_win = Tk()
    fit_win.title("Import and fit data")

    # Return chi^2 for given smeared spec
    def chi_square(smear_spec):
        # Empty list of energies matching imported spec
        energies_series = pd.Series(
            np.nan, index=import_dat.index
        )
        # Don't double count the smeared spec
        energies_series = energies_series[
            ~energies_series.index.isin(smear_spec.index)]
        # Concat with smeared spec
        inter_smear_dat = pd.concat([smear_spec, energies_series])
        inter_smear_dat.sort_index(inplace=True)
        # Interpolate to fill the new values between and AFTER import points
        inter_smear_dat.interpolate(method="linear", 
            limit_direction="both", 
            inplace=True)

        # Get rid of extra points we don't want from original smear spec
        inter_smear_dat = inter_smear_dat[
            inter_smear_dat.index.isin(import_dat.index)]

        # Normalise smeared plot
        inter_smear_dat_int = np.trapz(inter_smear_dat)
        inter_smear_dat_norm = inter_smear_dat.apply(
            lambda x: x/inter_smear_dat_int)
        diff_dat = import_dat_norm.subtract(inter_smear_dat_norm)
        diff_sq_dat = diff_dat.apply(lambda x: x**2)

        # print(diff_sq_dat.sum())
        # inter_smear_dat_norm.plot()
        # import_dat_norm.plot()
        # plt.show()

        return diff_sq_dat.sum()

    def fit_data(*args):
        def fit_recursive(fit_check_var_index=0):
            # If there are no more parameters to fit
            if(fit_check_var_index >= len(fit_check_vars)):
                # Calc chi_square, append parameters to df and return
                return
            # If this parameter needs to be fit
            elif(fit_check_vars[fit_check_var_index].get()):
                for i in range(N_STEPS):
                    # SET THIS PARAM TO NEW, PREDEFINED VALUE
                    fit_recursive(fit_check_var_index+1)
            else:
                # Leave this parameter as it is, move onto next
                fit_recursive(fit_check_var_index+1)
            return

        total_int_spec = pd.Series(0, index=ENERGIES)

        # Start at minimum values of defined range, calc chi-square
        for reactor in reactors:
            osc_spec = reactor.osc_spec(
                dm_21=(DM_21_FIT-DM_21_RANGE), 
                s_2_12=math.sin(2*(THET_12_FIT-THET_12_RANGE))**2, 
                dm_31=(DM_31_FIT-DM_31_RANGE),
                s_2_13=math.sin(2*(THET_13_FIT-THET_13_RANGE))**2, 
                period=period
            )
            int_spec = reactor.int_spec(osc_spec, "nu")
            total_int_spec = total_int_spec.add(int_spec)
        smear_spec = wit_smear.smear(total_int_spec, "nu")
        # Takes e+ spec as input so need to offset smear to match
        smear_spec = smear_spec.rename(UNDO_OFFSET_UP_DICT)

        # List of parameter values and their chi-squared to the data
        # DON'T START PARAMS THAT WONT BE FIT AT MIN
        fit_dat = [
            [DM_21_FIT-DM_31_RANGE,
                THET_12_FIT-THET_12_RANGE,
                DM_31_FIT-DM_31_RANGE,
                THET_13_FIT-THET_13_RANGE,
            chi_square(smear_spec)]]

        # Stores the param centre and range for a given cycle
        param_info = [
            [DM_21_FIT,THET_12_FIT,DM_31_FIT,THET_13_FIT],
            [DM_21_RANGE,THET_12_RANGE,DM_31_RANGE,THET_13_RANGE]
        ]
        
        for i in range(N_CYCLES):
            print("Fit cycle %i" % i)
            fit_recursive()
            # FILL LIST OF PARAMETER VALUES BASED ON CURRENT CYCLE
            # AN;D BEST FIT VALUE FOUND SO FAR IN THE DF
            best_fit_index = -1
            best_fit_chi = 1e6
            for i,row in enumerate(fit_dat):
                # Find index of minimum chi square so far
                if(row[-1] < best_fit_chi):
                    best_fit_index = i
            
            # Set new parameter centres to best fit params
            param_info[0] = param_info[i][:-1]
            # Shrink range to CYCLE_FACTOR of previous range
            param_info[1] = [x*CYCLE_FACTOR for x in param_info[1]] 

        return

    osc_fit_desc_label = Label(fit_win,
        text = "Choose which parameters to fit")
    osc_fit_desc_label.grid(column=0,row=2)

    # Selecting which osc vars to fit
    dm_21_fit_var = IntVar(value=0)
    dm_21_fit_check = Checkbutton(fit_win,
            text="delta m^2_21",
            variable=dm_21_fit_var)
    dm_21_fit_check.grid(column=0, row=3)

    thet_12_fit_var = IntVar(value=0)
    thet_12_fit_check = Checkbutton(fit_win,
            text="theta_12",
            variable=thet_12_fit_var)
    thet_12_fit_check.grid(column=0, row=4)

    dm_31_fit_var = IntVar(value=0)
    dm_31_fit_check = Checkbutton(fit_win,
            text="delta m^2_31",
            variable=dm_31_fit_var)
    dm_31_fit_check.grid(column=0, row=5)

    thet_13_fit_var = IntVar(value=0)
    thet_13_fit_check = Checkbutton(fit_win,
            text="theta_13",
            variable=thet_13_fit_var)
    thet_13_fit_check.grid(column=0, row=6)


    fit_check_vars = [
        dm_21_fit_var,
        thet_12_fit_var,
        thet_13_fit_var,
        dm_31_fit_var
    ]

    fit_button = Button(fit_win, text="Fit", command=fit_data)
    fit_button.grid(column=0, row=7)

