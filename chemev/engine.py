"""
engine
======

The file that contains the functions that run the chemical evolution model.
"""

import numpy as np, pickle
import logging
from . import zones

simfile = 'g1536_sfh+gas.pck'
def run(n_disk_zones=1,start_time=0,end_time=13.73e9,time_step=1e7,
        sf_mode='sim'):

    disk_zones = zones.create_zones(n_disk_zones)

    time_steps = np.arange(start_time,end_time,time_step)

    if sf_mode == 'sim':
    #Read in SFH and disk gas mass evolution
        sfh_dat = pickle.load(open(simfile))
        sfh = np.interp(time_steps,sfh_dat['startimes'],sfh_dat['sfh'],left=0)
        gas_mass_ev = np.interp(time_steps,sfh_dat['gastimes'],sfh_dat['gasmass'],left=0)
    else: 
        sf_mode = 'kennicutt'

    # Read in chemical evolution tables
    chemev_dat = pickle.load(open('time_yields.pck'))


    for it,time in enumerate(time_steps):
        print 'Time: %g'%time

        gas_mass = gas_mass_ev[it]
        sf_mass = sfh[it]
        for zone in disk_zones:
            zone.enrich(time)
            zone.form_star(time,sf_mode=sf_mode,mass=sf_mass)
            
