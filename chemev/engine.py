"""
engine
======

The file that contains the functions that run the chemical evolution model.
"""

import numpy as np, pickle, pyfits, os
try:
    import pynbody
    have_pynbody = True
except:
    have_pynbody = False
import logging, pdb
from . import zones, enrich

simfile = 'g1536_sfh+gas.pck'
def run(n_disk_zones=1,start_time=0,end_time=13.73e9,time_step=1e7,
        max_disk_r=25,disk_gas_mass=1e10,h_r=4,
        sf_mode='sim',init_abunds={'O':0,'Fe':0,'Mg':0},
        snii_yields='kobayashi',snia_yields='iwamoto',agb_yields='karakas',
        outfits='simstar.fits'):

    star_names = ['tform','init_mass','mass','Z','zone']
    star_names.extend(init_abunds.keys())

    if not have_pynbody:
    # Dynamically create star type depending on how many elemental abundances
    # are being tracked.
        st_formats = ['f','f','f','f']
        st_formats.extend(len(init_abunds)*['f'])
        zones.star_type = np.dtype({'names':star_names,'formats':st_formats})
        
    # Set up geometry
    max_disk_r, h_r = float(max_disk_r), float(h_r)
    disk_zones = zones.create_disk_zones(n_disk_zones,init_abunds=init_abunds,
                                         max_disk_r=max_disk_r,
                                         disk_gas_mass=disk_gas_mass,h_r=h_r)

    # Set up timesteps
    time_steps = np.arange(start_time,end_time,time_step)
    n_t_steps = len(time_steps)

    stars = pynbody.snapshot.SimSnap()
    stars._num_particles = n_disk_zones*n_t_steps
    stars._filename = "<created>"
    stars._create_arrays(star_names,1)
    stars._family_slice['star'] = slice(0,n_disk_zones*n_t_steps)
    for iz,zone in enumerate(disk_zones):
        stars[iz*n_t_steps:(iz+1)*n_t_steps]['zone'] = iz
        stars[iz*n_t_steps:(iz+1)*n_t_steps]['tform'] = time_steps

    # Set up AGB and SNII enrichment models
    mydir=os.path.dirname(__file__)
    enrich.make_table('agb',time_steps,
                      infile=mydir+'/yields/agb/'+agb_yields+'.pck')
    enrich.make_table('snii',time_steps,
                      infile=mydir+'/yields/snii/'+snii_yields+'.pck')
    enrich.make_snia_table(time_steps,
                           infile=mydir+'/yields/snia/'+snia_yields+'.pck')

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
        sf_mass_tot = sfh[it]
        Sigma0 = gas_mass / (np.pi*2*h_r*(h_r - 
                                np.exp(-max_disk_r/h_r)*
                                (h_r+max_disk_r)))
        prefac = np.pi*2.0*h_r*Sigma0
        
        # All zones will have the same star ages as the first one.
        ages = time - stars[0:it]['tform']

        for iz,zone in enumerate(disk_zones):
            izstar = np.arange(iz*n_t_steps,iz*n_t_steps+it)

            zone_gas_mass = prefac*(np.exp(-zone.min_r/h_r)*(h_r+zone.min_r) - 
                                np.exp(-zone.max_r/h_r)*(h_r+zone.max_r))

            rel_masses = zone.enrich(ages,stars[izstar]['Z'],
                                     stars[izstar]['mass'],
                                     time_step,zone_gas_mass)
            stars[izstar]['mass']-=stars[izstar]['mass']*rel_masses*time_step

            sf_mass = sf_mass_tot * zone.mass / gas_mass
            if sf_mass >0: 
                # Initialize new star
                istar = iz*n_t_steps+it
                stars[istar]['init_mass'] = sf_mass
                stars[istar]['mass'] = sf_mass
                stars[istar]['Z'] = zone.Z
                for el in init_abunds:
                    stars[istar][el] = zone.abunds[el]

    stars.write(filename='stars.tips',fmt=pynbody.snapshot.tipsy.TipsySnap)

    return stars,disk_zones
