#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

""" 
    The portion of SKReact dealing with neutrino production in
    nuclear reactors.
"""

import params
import pandas
from math import sin, cos, tan, sqrt, radians
from calendar import monthrange

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
            lf_monthly):

        self.country = country
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.core_type = core_type
        self.mox = mox
        self.p_th = p_th # MW
        self.lf_monthly = lf_monthly #Pandas series

    def n_nu(self, period = "Max"):
        nu_per_s = self.p_th*params.NU_PER_MW 
        if(period == "Max" or period == "max"): # Yearly at reference P
            return 365*24*60*nu_per_s
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
                    lf_month = self.lf_monthly[
                            "LF_"+str(year)+"/"+str(month).zfill(2)]

                    n_nu_month = (n_days_in_month*24*60)
                    n_nu_month *= (lf_month*self.p_th)
                    n_nu_month *= params.NU_PER_MW

                    n_nu_tot += n_nu_month

            return n_nu_tot
        elif(len(period) == 7): # Specific month YYYY/MM
            year  = int(period[:4])
            month = int(period[5:])
            n_days_in_month = monthrange(year,month)[1]
            lf_month = self.lf_monthly["LF_"+str(year)+"/"+str(month).zfill(2)]
            n_nu_month = (n_days_in_month*24*60)
            n_nu_month *= (lf_month*self.p_th)
            n_nu_month *= params.NU_PER_MW
            return n_nu_month
        elif(len(period) == 4): # Specific year YYYY
            year = int(period)

            # Cycle through all months calculating nu per month
            n_nu_tot = 0
            for month in range(1,13):
                n_days_in_month = monthrange(year,month)[1]
                # Query the specific month from the LF series
                lf_month = self.lf_monthly[
                        "LF_"+str(year)+"/"+str(month).zfill(2)]

                n_nu_month = (n_days_in_month*24*60)
                n_nu_month *= (lf_month*self.p_th)
                n_nu_month *= params.NU_PER_MW

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
        a = params.EARTH_R_EQUATOR**2*cos(latitude)
        b = params.EARTH_R_POLAR**2*sin(latitude)
        c = params.EARTH_R_EQUATOR*cos(latitude)
        d = params.EARTH_R_POLAR*sin(latitude)

        r = sqrt((a*a + b*b)/(c*c + d*d))

        return r

    """
        Returns sin of geocentric latitude from geodetic latitude
    """
    def sin_geocentric(self, latitude):
        tan_a = params.EARTH_R_POLAR*tan(latitude)/params.EARTH_R_EQUATOR
        sin_a = tan_a/sqrt(1+tan_a*tan_a)
        return sin_a

    """
        Returns cos of geocentric latitude from geodetic latitude
    """
    def cos_geocentric(self, latitude):
        tan_a = params.EARTH_R_POLAR*tan(latitude)/params.EARTH_R_EQUATOR
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

        lat_sk_rad = radians(params.SK_LAT)
        long_sk_rad = radians(params.SK_LONG)

        r_react = self.dist_to_earth_centre(lat_react_rad)
        r_sk    = self.dist_to_earth_centre(lat_sk_rad) + params.SK_ALT

        x_react = r_react*self.cos_geocentric(lat_react_rad)*cos(long_react_rad)
        y_react = r_react*self.cos_geocentric(lat_react_rad)*sin(long_react_rad)
        z_react = r_react*self.sin_geocentric(lat_react_rad)

        x_sk = r_sk*self.cos_geocentric(lat_sk_rad)*cos(long_sk_rad)
        y_sk = r_sk*self.cos_geocentric(lat_sk_rad)*sin(long_sk_rad)
        z_sk = r_sk*self.sin_geocentric(lat_sk_rad)

        dist = sqrt((x_react-x_sk)**2 + (y_react-y_sk)**2 + (z_react-z_sk)**2)
        
        return dist
