#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alex Goldsack"

""" 
    The portion of SKReact dealing with neutrino production in
    nuclear reactors.
"""

import params
from math import sin, cos, tan, sqrt, radians

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
            p_monthly):

        self.country = country
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.core_type = core_type
        self.mox = mox
        self.p_th = p_th
        self.p_monthly = p_monthly.copy()

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
