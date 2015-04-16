"""
engine
======

The file that contains the functions that run the chemical evolution model.
"""

import numpy as np, pickle, pdb
import logging
from . import zones, data

simfile = 'g1536_sfh+gas.pck'
def run(n_disk_zones=1,start_time=0,end_time=13.73e9,time_step=1e7,
        sf_mode='sim',init_abunds={'O':0,'Fe':0,'Mg':0}):

    data.star_type = np.dtype({'names':['tform','init_mass','mass',
                                        'Z'].extend(init_abunds.keys()),
                               'formats':['f','f','f',
                                          'f'].extend(len(init_abunds)*'f')})
    disk_zones = zones.create_disk_zones(n_disk_zones,init_abunds=init_abunds)

    time_steps = np.arange(start_time,end_time,time_step)

    if sf_mode == 'sim':
    #Read in SFH and disk gas mass evolution
        sfh_dat = pickle.load(open(simfile))
        sfh = np.zeros(len(time_steps))
        for it,t_step in enumerate(time_steps[:-1]):
            in_range = (sfh_dat['startimes']*1e9 > t_step) & \
                (sfh_dat['startimes']*1e9 < time_steps[it+1])
            sfh[it] = sfh_dat['sfh'][in_range].sum()
        gas_mass_ev = np.interp(time_steps,
                                np.array(sfh_dat['gastimes'])*1e9,
                                np.array(sfh_dat['gasmass']),left=0)
    else: 
        sf_mode = 'kennicutt'

    for it,time in enumerate(time_steps):
        print 'Time: %g'%time

        gas_mass = gas_mass_ev[it]
        sf_mass = sfh[it]

        for zone in disk_zones:
            zone.enrich(time,time_step,gas_mass)
            if sf_mass >0: zone.form_star(time,sf_mode=sf_mode,mass=sf_mass)
            
    import pdb; pdb.set_trace()
