"""
engine
======

The file that contains the functions that run the chemical evolution model.
"""

import numpy as np, pickle
import logging
from . import zones

def run(n_disk_zones=1,start_time=0,end_time=13.73e9,time_step=1e6):

    #Read in SFH and disk gas mass evolution
    sfh_dat = pickle.load(open('g1536_sfh+gas.pck'))

    # Read in chemical evolution tables
    chemev_dat = pickle.load(open('time_yields.pck'))

    disk_zones = zones.create_zones(n_disk_zones)

    time_steps = np.arange(start_time,end_time,time_step)

    for time in time_steps:
        print 'Time: %g'%time
        for zone in disk_zones:
            zone.form_star(time)
            zone.enrich()
            
