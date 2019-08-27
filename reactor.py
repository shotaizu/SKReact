#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

""" 
The portion of SKReact dealing with neutrino production in
nuclear reactors and subsequent oscillation.
"""

import params
from params import *
import pandas
from math import sin, cos, tan, sqrt, radians
from calendar import monthrange
import numpy as np
import math

class Reactor:

    # Initialiser
    def __init__(self,
            country,
            name,
            latitude,
            longitude,
            core_type,
            mox,
            p_th,
            lf_monthly,
            default=True):

        self.country = country
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        #core_type is checked later, need to remove whitespace
        self.core_type = core_type.rstrip()
        self.mox = mox
        self.p_th = p_th # MW
        self.dist_to_sk = self.dist_to_sk()
        self.lf_monthly = lf_monthly #Pandas series
        self.p_monthly = self.p_monthly()
        self.p_r_monthly = self.p_r_monthly()
        self.default = default # If the reactor came from the xls

    # Monthly power output calculate from load factor and p_th
    def p_monthly(self):
        # Same format as lf
        index = self.lf_monthly.index 
        p_list = [self.pth*lf for lf in self.lf_monthly.tolist()]
        return pd.Series(p_list, index=index)

    # Monthly power/r^2 output calculate from p_monthly and dist_to_sk
    def p_r_monthly(self):
        # Same format as lf
        index = self.p_monthly.index 
        p_r_list = [p/(self.dist_to_sk**2) for p in self.p_monthly.tolist()]
        return pd.Series(p_r_list, index=index)

    def add_to_lf(self, date, lf):
        # print(self.name)
        # print(date)
        # print(lf)
        self.lf_monthly.set_value(date, lf)
        self.p_monthly.set_value(date, lf*self.p_th)
        self.p_r_monthly.set_value(date, lf*self.p_th/(self.dist_to_sk**2))
        return

    # Calculate the number of neutrinos produced in given period
    # TODO: Move the common calcs outside the if statement
    def n_nu(self, period = "Max"):
        # Pre-calculating the nu per second at reference power for self
        nu_per_s = self.p_th*NU_PER_MW 
        if(period == "Max" or period == "max"): # Yearly at reference P
            return 365*24*60*60*nu_per_s
        elif(len(period) == 15): # Inclusive period YYYY/MM-YYYY/MM
            year_start  = int(period[:4])
            month_start = int(period[5:7])
            year_end  = int(period[8:12])
            month_end = int(period[13:])

            # Cycle through all months calculating nu per month
            month_range_start = month_start
            month_range_end = 13
            n_nu_tot = 0
            for year in range(year_start,year_end+1):
                # Start from Jan after first year
                if(year != year_start):
                    month_range_start = 1
                # Only go up to end of period in final year
                if(year == year_end):
                    month_range_end = month_end+1 # For inclusivity
                for month in range(month_range_start,month_range_end):
                    n_days_in_month = monthrange(year,month)[1]
                    # Query the specific month from the LF series
                    lf_month = self.lf_monthly["%i/%02i" % (year, month)]
                    lf_month /= 100 #To be a factor, not %age
                    n_nu_month = (n_days_in_month*24*60*60)
                    n_nu_month *= (lf_month*nu_per_s)

                    n_nu_tot += n_nu_month

            return n_nu_tot
        elif(len(period) == 7): # Specific month YYYY/MM
            year  = int(period[:4])
            month = int(period[5:])
            n_days_in_month = monthrange(year,month)[1]
            lf_month = self.lf_monthly["%i/%02i" % (year, month)]
            lf_month /= 100
            n_nu_month = (n_days_in_month*24*60*60)
            n_nu_month *= (lf_month*nu_per_s)
            return n_nu_month
        elif(len(period) == 4): # Specific year YYYY
            year = int(period)

            # Cycle through all months calculating nu per month
            n_nu_tot = 0
            for month in range(1,13):
                n_days_in_month = monthrange(year,month)[1]
                # Query the specific month from the LF series
                lf_month = self.lf_monthly["%i/%02i" % (year, month)]
                lf_month /= 100

                n_nu_month = (n_days_in_month*24*60*60)
                n_nu_month *= (lf_month*n_nu_s)

                n_nu_tot += n_nu_month

            return n_nu_tot
        else:
            print("reactor.no_nu() requires either YYYY/MM, YYYY or \"Max\" "
                "(per year) for period of nu production.")
            exit()

    """ 
    Earth bulges a the equator, this gives distance to
    centre of the Earth as a function of latitude
    """
    def dist_to_earth_centre(self, latitude):
        a = EARTH_R_EQUATOR**2*cos(latitude)
        b = EARTH_R_POLAR**2*sin(latitude)
        c = EARTH_R_EQUATOR*cos(latitude)
        d = EARTH_R_POLAR*sin(latitude)

        r = sqrt((a*a + b*b)/(c*c + d*d))

        return r

    """
    Returns sin of geocentric latitude from geodetic latitude
    """
    def sin_geocentric(self, latitude):
        tan_a = EARTH_R_POLAR*tan(latitude)/EARTH_R_EQUATOR
        sin_a = tan_a/sqrt(1+tan_a*tan_a)
        return sin_a

    """
    Returns cos of geocentric latitude from geodetic latitude
    """
    def cos_geocentric(self, latitude):
        tan_a = EARTH_R_POLAR*tan(latitude)/EARTH_R_EQUATOR
        cos_a = 1/sqrt(1+tan_a*tan_a)
        return cos_a

    """
    Use Lat and Long info to calc distance to SK in km
    Assume reactors are at sea level, very reasonable
    assumption, given most are on coastline
    """
    def dist_to_sk(self):
        lat_react_rad = radians(self.latitude)
        long_react_rad = radians(self.longitude)

        lat_sk_rad = radians(SK_LAT)
        long_sk_rad = radians(SK_LONG)

        r_react = self.dist_to_earth_centre(lat_react_rad)
        r_sk    = self.dist_to_earth_centre(lat_sk_rad) + SK_ALT

        x_react = r_react*self.cos_geocentric(lat_react_rad)*cos(long_react_rad)
        y_react = r_react*self.cos_geocentric(lat_react_rad)*sin(long_react_rad)
        z_react = r_react*self.sin_geocentric(lat_react_rad)

        x_sk = r_sk*self.cos_geocentric(lat_sk_rad)*cos(long_sk_rad)
        y_sk = r_sk*self.cos_geocentric(lat_sk_rad)*sin(long_sk_rad)
        z_sk = r_sk*self.sin_geocentric(lat_sk_rad)

        dist = sqrt((x_react-x_sk)**2 + (y_react-y_sk)**2 + (z_react-z_sk)**2)
        
        return dist

    """
    Getting E spectrum from 5th order polynomial for given coefficients
    """
    def f_from_poly(self,energy,coeffs):
        flux = 0
        for i,a in enumerate(coeffs):
            flux += a*energy**i
        return math.exp(flux)

    """
    Use 5th order polynomial to simulate E spectrum produced by reactors,
    taking into account reactor type and whether the reactor uses MOX
    """
    def e_spectra(self):
        energies = np.linspace(E_MIN, E_MAX, E_BINS)
        core_type = self.core_type
        if(self.mox):
            core_type = "MOX"

        # Fuel fractions for this type of core
        u_235_frac = FUEL_MAKEUP.loc[core_type]["U_235"]
        pu_239_frac = FUEL_MAKEUP.loc[core_type]["Pu_239"]
        u_238_frac = FUEL_MAKEUP.loc[core_type]["U_238"]
        pu_241_frac = FUEL_MAKEUP.loc[core_type]["Pu_241"]

        u_235_spectrum = [(u_235_frac/U_235_Q)*self.f_from_poly(energy,U_235_A)
                for energy in energies]
        pu_239_spectrum = [(pu_239_frac/PU_239_Q)*self.f_from_poly(energy,PU_239_A) 
                for energy in energies]
        u_238_spectrum = [(u_238_frac/U_238_Q)*self.f_from_poly(energy,U_238_A) 
                for energy in energies]
        pu_241_spectrum = [(pu_241_frac/PU_241_Q)*self.f_from_poly(energy,PU_241_A) 
                for energy in energies]
        tot_spectrum = [sum(f) for f in zip(u_235_spectrum, 
            pu_239_spectrum, 
            u_238_spectrum, 
            pu_241_spectrum)]

        spectrum_data = { "U_235": u_235_spectrum,
                "Pu_239": pu_239_spectrum,
                "U_238": u_238_spectrum,
                "Pu_241": pu_241_spectrum,
                "Total" : tot_spectrum}

        e_spectra = pd.DataFrame(spectrum_data, index=energies)

        return e_spectra

    """
    Calculating the spectrum of ALL oscillated nu E at SK
    TODO: Add in hierarchy support (I think it barely changes it)
    """
    def oscillated_spec(self,
            dm_21 = DM_21,
            c_13 = C_13_NH,
            s_2_12 = S_2_12,
            s_13 = S_13_NH,
            period = "Max"):

        # Finding total load factor 
        year_start  = int(period[:4])
        month_start = int(period[5:7])
        year_end  = int(period[8:12])
        month_end = int(period[13:])

        # Cycle through all months summing load factor*t 
        lf_s_sum = 0
        month_range_start = month_start
        month_range_end = 13
        n_nu_tot = 0
        for year in range(year_start,year_end+1):
            # Start from Jan after first year
            if(year != year_start):
                month_range_start = 1
            # Only go up to end of period in final year
            if(year == year_end):
                month_range_end = month_end+1 # For inclusivity
            for month in range(month_range_start,month_range_end):
                n_days_in_month = monthrange(year,month)[1]
                # Query the specific month from the LF series
                # print(self.lf_monthly)
                try:
                    lf_month = float(self.lf_monthly["%i/%02i" % (year, month)])
                except TypeError:
                    print("Load factor data for reactor "
                            + self.name
                            + " in month %i/%02i" % (year,month)
                            + " not float compatible")
                    exit()
                except KeyError:
                    print("Error with " 
                            + self.name 
                            + " in or around file DB%i.xls" % year)
                    print("Does not have entry for this year.")
                    exit()
                lf_month /= 100 #To be a factor, not %age
                lf_s_sum += lf_month*n_days_in_month*24*60*60

        spec_pre_factor = self.p_th*lf_s_sum
        # From PHYSICAL REVIEW D 91, 065002 (2015)
        # E in MeV, l in km
        p_ee = lambda e: c_13*c_13*(1-s_2_12*(math.sin(1.27*dm_21*l*1e3/e))**2)+s_13*s_13

        l = self.dist_to_sk
        energies = np.linspace(E_MIN, E_MAX, E_BINS)
        # Don't think I'll need osc. spectra of individual fuels
        e_spec = self.e_spectra()["Total"].tolist()
        # Use max to skip 0
        # osc_e_spec = [f*p_ee(max(e,1e-5)) for f,e in zip(e_spec,energies)]
        osc_e_spec = []
        for f,e in zip(e_spec,energies):
            if(e > IBD_MIN):
                # Can't forget to drop off with r^2
                osc_e_spec.append(spec_pre_factor*f*p_ee(e)/(self.dist_to_sk**2))
            else:
                osc_e_spec.append(0)
        osc_spec = pd.Series(osc_e_spec, index=energies)
        return osc_spec

    """
    Spectrum of INCIDENT oscillated nu E at SK
    Takes oscillated spec and multiplies by xsec
    """
    def incident_spec(self,
            dm_21 = DM_21,
            c_13 = C_13_NH,
            s_2_12 = S_2_12,
            s_13 = S_13_NH,
            period = "Max"):
        # From PHYSICAL REVIEW D 91, 065002 (2015)
        e_e = lambda e: e - DEL_NP
        p_e = lambda e: math.sqrt(e_e(e)**2 - M_E*M_E)
        e_exp = lambda e: e**(-0.07056+0.02018*math.log(e)-0.001953*(math.log(e))**3)
        xsec = lambda e: 1e-43*p_e(e)*e_e(e)*e_exp(e) # cm^2

        energies = np.linspace(E_MIN, E_MAX, E_BINS)
        osc_spec = self.oscillated_spec(dm_21,c_13,s_2_12,s_13,period)
        incident_spec_dat = []
        for e,f in zip(energies,osc_spec):
            if(e > IBD_MIN):
                incident_spec_dat.append(f*xsec(e))
            else:
                incident_spec_dat.append(0)

        incident_spec = pd.Series(incident_spec_dat, index=energies)

        return incident_spec
