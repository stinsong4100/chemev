"""
FILE
    iwamoto99_yields.py

DESCRIPTION

    Convert Iwamoto et al. (1999) W70 SNIa yields ASCII data table to pickled
    arrays.
    
"""

import numpy as np
import os
import pickle
import string


def pickle_write(obj, filename):
    '''write pickled object'''
    fout = open(filename, 'w')
    pickle.dump(obj, fout, 2)
    fout.close()


version = 0.0

##### Load Data #####
home = os.path.expanduser('~')
stem_flexce = home + '/FlexCE_v' + str(version) + '/'
stem_yields = stem_flexce + 'yields/'
stem_yldgen = stem_yields + 'general/'
stem_i99 = stem_yields + 'iwamoto99/'
####################


### Iwamoto et al. (1999) W70 SNIa yields ###

'''models:
w7: single-degenerate model from Nomoto et al. (1984)
w70: zero-metallicity version of w7
wdd1: 
wdd2:
wdd3: 
cdd1:
cdd2:
'''

# Read in yields
snia_in = stem_i99 + 'iwamoto99.txt'
w7, w70, wdd1, wdd2, wdd3, cdd1, cdd2 = \
     np.loadtxt(snia_in, usecols=(2,3,4,5,6,7,8), unpack=True)
snia_dict = {'w7': w7, 'w70': w70, 'wdd1': wdd1, 'wdd2': wdd2, 'wdd3': wdd3,
             'cdd1': cdd1, 'cdd2': cdd2}
snia_sym_in = np.loadtxt(snia_in, usecols=(0,), dtype=str, unpack=True)
snia_sym = np.array([''.join((item[2:], item[:2])) for item in snia_sym_in])

# Read in isotopes
species = np.loadtxt(stem_yldgen + 'species.txt', dtype=str, usecols=(1,))
n_species = len(species)
el_name = np.array([item.strip(string.digits) for item in species])
iso_mass = np.array([int(item.strip(string.ascii_letters))
                     for item in species])
elements = []
for item in el_name:
    if item not in elements:
        elements.append(item)

elements = np.array(elements)


# Solar abundances are prevalence "by mass"
solar_iso = np.loadtxt(stem_yldgen + 'Solar_isotopes.txt', skiprows=1,
                       dtype=str, usecols=(0,))
solar_ab = np.loadtxt(stem_yldgen + 'Solar_isotopes.txt', skiprows=1,
                      usecols=(1,))
solar_el_name = np.array([item.strip(string.digits) for item in solar_iso])


# indices within "species" array of the elements for which CL04 give a solar
# abundance
ind_iso = []
for i in xrange(len(solar_iso)):
    ind_iso.append(np.where(solar_iso[i] == species)[0][0])

ind_iso = np.array(ind_iso)


# map elemental yields onto the dominant isotope
snia_yields = {}
for k in snia_dict.iterkeys():
    snia_yields[k] = np.zeros(n_species)
    for j in xrange(n_species):
        if species[j] in snia_sym:
            snia_yields[k][j] = snia_dict[k][np.where(snia_sym == species[j])]


# write to file
for k in snia_yields.iterkeys():
    pickle_write(snia_yields[k], stem_i99 + k + '_yields.pck')
