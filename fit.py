from params import *
from reactor import Reactor
from smear import Smear
from tkinter import *
from tkinter import filedialog
import tkinter.ttk as ttk
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import math
import os
import sys


def main():
    if len(sys.argv) < 2:
        print("NO FILE PREFIX TO PULL FROM")
        print("Either run through SKReact or give file or file prefix to fit")
        return

    import_filename = sys.argv[1]
    try:
        period = sys.argv[2]
    except:
        period = None

    try:
        out_filename = sys.argv[3]
    except:
        out_filename = None

    # Get the reactors
    with open(REACT_PICKLE, "rb") as pickle_file:
        reactors = pickle.load(pickle_file)

    # Set their spectra for current skreact params
    for reactor in reactors:
        reactor.set_all_spec()

    # And the smearing matrix
    try:
        wit_smear = Smear(WIT_SMEAR_FILE)
    except FileNotFoundError:
        print("Smear file " + WIT_SMEAR_FILE + " not found!")
        print("Cannot import smearing information.")
        return

    fit_win(import_filename, reactors, period, wit_smear, out_filename)


# Window to handle data import and fitting
def fit_win(import_filename, reactors, period, wit_smear, out_filename=None):

    # Now make window
    # fit_win = Toplevel(skreact_win)
    fit_win = Tk()
    fit_win.title("Import and fit data")

    # Return chi^2 for given smeared spec
    def chi_square(smear_spec, import_dat, plot_norm=False):
        # Empty list of energies matching imported spec
        energies_series = pd.Series(np.nan, index=import_dat.index)
        # Don't double count the smeared spec
        energies_series = energies_series[~energies_series.index.isin(smear_spec.index)]
        # Concat with smeared spec
        inter_smear_dat = pd.concat([smear_spec, energies_series])
        inter_smear_dat.sort_index(inplace=True)
        # Interpolate to fill the new values between and AFTER import points
        inter_smear_dat.interpolate(limit_direction="both", inplace=True)

        # Get rid of extra points we don't want from original smear spec
        inter_smear_dat = inter_smear_dat[inter_smear_dat.index.isin(import_dat.index)]

        # Normalise smeared plot
        inter_smear_dat_int = np.trapz(inter_smear_dat)
        inter_smear_dat_norm = inter_smear_dat.apply(lambda x: x / inter_smear_dat_int)
        import_dat_int = np.trapz(import_dat)
        import_dat_norm = import_dat.apply(lambda x: x / import_dat_int)
        import_dat_norm_int = np.trapz(import_dat_norm)
        diff_dat = import_dat_norm.subtract(inter_smear_dat_norm)
        diff_sq_dat = diff_dat.apply(lambda x: x ** 2)

        if plot_norm:
            print(diff_sq_dat.sum())
            # import_dat_norm.plot.bar(width=1.0)
            import_dat_norm.plot()
            inter_smear_dat_norm.plot()
            # plt.legend()
            plt.show()

        return diff_sq_dat.sum()

    def fit_data(import_filename):
        try:
            import_dat_df = pd.read_csv(import_filename)
        except:
            print("Cannot read import data!")
            return

        # Change to series
        import_dat_df.columns = ["energy", "bin_content"]
        import_dat = import_dat_df.set_index("energy")["bin_content"]

        # Area normalise
        import_dat_int = np.trapz(import_dat)
        import_dat_norm = import_dat.apply(lambda x: x / import_dat_int)

        import_dat_norm_int = np.trapz(import_dat_norm)

        # Calculate the smeared spec for the current parameters
        def calc_smear():
            total_int_spec = pd.Series(0, index=ENERGIES)
            for reactor in reactors:
                osc_spec = reactor.osc_spec(
                    dm_21=param_values[0],
                    s_2_12=math.sin(2 * param_values[1]) ** 2,
                    dm_31=param_values[2],
                    s_2_13=math.sin(2 * param_values[3]) ** 2,
                    period=period,
                )
                int_spec = reactor.int_spec(osc_spec, "nu")
                total_int_spec = total_int_spec.add(int_spec)
            smear_spec = wit_smear.smear(total_int_spec, "nu")
            # Takes e+ spec as input so need to offset smear to match
            return smear_spec.rename(UNDO_OFFSET_UP_DICT)

        # Cycle through all parameters space
        def fit_recursive(param_index=0):
            nonlocal fit_check_vars
            # If there are no more parameters to fit
            if param_index >= len(fit_check_vars):
                # Calc chi_square, append parameters to df and return
                fit_dat.append(param_values + [chi_square(calc_smear(), import_dat)])
                print(fit_dat[-1])
                # plt.cla()
                # plt.plot([row[0] for row in fit_dat[prev_cycle_n_rows:]],
                #     [row[-1] for row in fit_dat[prev_cycle_n_rows:]])
                # plt.pause(0.05)
                return
            # If this parameter needs to be fit
            elif fit_check_vars[param_index].get():
                this_best_fit_chi = 1e6
                print("Cycling param %i" % param_index)
                # for i in range(n_steps):
                step = 0
                while fit_dat[-1][-1] <= this_best_fit_chi or step < N_STEPS:
                    this_best_fit_chi = fit_dat[-1][-1]
                    best_fit_index = len(fit_dat) - 1
                    fit_recursive(param_index + 1)
                    # Step the param value forward
                    param_values[param_index] += param_cycle_info[2][param_index]
                    step += 1
            else:
                # Leave this parameter as it is, move onto next
                fit_recursive(param_index + 1)
            return

        # Stores the param centre and range for a given cycle
        param_cycle_info = [
            [DM_21_FIT, THET_12_FIT, DM_31_FIT, THET_13_FIT],
            [DM_21_RANGE, THET_12_RANGE, DM_31_RANGE, THET_13_RANGE],
            [
                2 * DM_21_RANGE / N_STEPS,
                2 * THET_12_RANGE / N_STEPS,
                2 * DM_31_RANGE / N_STEPS,
                2 * THET_13_RANGE / N_STEPS,
            ],
        ]
        # List of param values to calc spec for at any one time
        param_values = [DM_21_FIT, THET_12, DM_31, THET_13]
        # Set min values

        for i, (fit, rng) in enumerate(zip(param_cycle_info[0], param_cycle_info[1])):
            param_values[i] = fit
            # If this one is to be fit, set it to min in range
            if fit_check_vars[i].get():
                param_values[i] -= rng

        # List of parameter values and their chi-squared to the data
        fit_dat = [param_values + [chi_square(calc_smear(), import_dat)]]

        best_fit_index = -1
        best_fit_chi = 1e6
        prev_cycle_n_rows = 0
        for i in range(N_CYCLES):
            print("Fit cycle %i" % i)
            fit_recursive()
            # Find index of minimum chi square so far
            # Skipping previous cycles
            best_fit_index = -1
            for i, row in enumerate(fit_dat[prev_cycle_n_rows:]):
                print(row)
                if row[-1] < best_fit_chi:
                    # Update best fit info
                    best_fit_index = i
                    best_fit_chi = row[-1]
            best_fit_index += prev_cycle_n_rows
            print("New best fit with chi-square %f" % best_fit_chi)
            print("With values:")
            print(fit_dat[best_fit_index][:-1])
            # Set new parameter centres to best fit params
            param_cycle_info[0] = fit_dat[best_fit_index][:-1]
            # Shrink range to CYCLE_FACTOR of previous range
            param_cycle_info[1] = [x * CYCLE_FACTOR for x in param_cycle_info[1]]
            # Update step size, double to go from -ve to +ve range
            param_cycle_info[2] = [2 * x / N_STEPS for x in param_cycle_info[1]]
            prev_cycle_n_rows = len(fit_dat)

            # Now set param values
            for i, (fit, rng) in enumerate(
                zip(param_cycle_info[0], param_cycle_info[1])
            ):
                param_values[i] = fit
                # If this one is to be fit, set it to min in range
                if fit_check_vars[i].get():
                    param_values[i] -= rng

        plt.close()
        print("Done fitting!")
        print("Final parameters:")
        print(fit_dat[best_fit_index][:-1])
        param_values = fit_dat[best_fit_index][:-1]
        chi_square(calc_smear(), import_dat, plot_norm=True)
        if out_filename != None:
            out_file = open(sys.argv[3], "a")
            out_file.write(",".join([str(x) for x in param_values]) + "\n")

        return

    osc_fit_desc_label = Label(fit_win, text="Choose which parameters to fit")
    osc_fit_desc_label.grid(column=0, row=2)

    # Selecting which osc vars to fit
    dm_21_fit_var = IntVar(value=1)
    dm_21_fit_check = Checkbutton(fit_win, text="delta m^2_21", variable=dm_21_fit_var)
    dm_21_fit_check.grid(column=0, row=3)

    thet_12_fit_var = IntVar(value=0)
    thet_12_fit_check = Checkbutton(fit_win, text="theta_12", variable=thet_12_fit_var)
    thet_12_fit_check.grid(column=0, row=4)

    dm_31_fit_var = IntVar(value=0)
    dm_31_fit_check = Checkbutton(fit_win, text="delta m^2_31", variable=dm_31_fit_var)
    dm_31_fit_check.grid(column=0, row=5)

    thet_13_fit_var = IntVar(value=0)
    thet_13_fit_check = Checkbutton(fit_win, text="theta_13", variable=thet_13_fit_var)
    thet_13_fit_check.grid(column=0, row=6)

    fit_check_vars = [dm_21_fit_var, thet_12_fit_var, thet_13_fit_var, dm_31_fit_var]

    # Decides to fit one or multiple files
    def fit_handler(*args):
        nonlocal period
        # Check if it's a prefix or actual file
        if import_filename[-4:] == ".csv":
            print("FITTING ONE FILE")
            if period == None:
                period = input("Enter period to fit for (YYYY/MM-YYYY/MM): ")
            fit_data(import_filename)
        else:
            print("FITTING MULTIPLE FILES")
            file_names = []
            for file in os.listdir(FIT_DIR):
                file_name = os.fsdecode(file)
                if file_name.startswith(import_filename):
                    file_names.append(file_name)

            file_names.sort()

            for file_name in file_names:
                print("Fitting file " + file_name)
                fit_data(file_name)

    fit_button = Button(fit_win, text="Fit", command=fit_handler)
    fit_button.grid(column=0, row=7)

    fit_win.mainloop()


if __name__ == "__main__":
    main()
