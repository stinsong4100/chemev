"""
FILE
    cescutti06_yields.py

DESCRIPTION

    Generates a finely spaced grid of Ba & Eu yields from SNII from the
    'empirical' r-process yields (determined by matching Ba and Eu abundances
    from a chemical evolution model that used these yields to observations) of
    Cescutti et al. (2001).

    Busso et al. (2001): M = 12--30 Msun; no metallicity dependence

"""

import numpy as np
import os
import pickle
from scipy import interpolate


def pickle_read(filename):
    '''read pickled object'''
    fin = open(filename, 'r')
    obj = pickle.load(fin)
    fin.close()
    return obj

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
stem_c06 = stem_yields + 'cescutti06/'
####################


##### Load Data ####
data_in = np.loadtxt(stem_c06 + 'cescutti06_rprocess.txt', unpack=True)
m_c06 = data_in[0]
x_c06 = data_in[1:]
c06_orig = x_c06 * m_c06 # convert X_Ba & X_Eu to a mass of Ba & Eu
####################


#####
# chemical evolution model mass bins
# IMF
alpha = 2.35
Gamma = 1. - alpha
alpha2 = 2. - alpha
m_min = 0.1
m_max = 100.
a = alpha2 / (m_max**alpha2 - m_min**alpha2)
m_cutoff = 8.

# Bins of Stars
'''Bin lower bounds (bins, bins_low, bins_high). Bin width (dbin_low,
dbin_high).  Number of bins (n_bins, n_bins_low, n_bins_high).  Average mass
per bin (m_ave_high, m_ave_low, m_ave), fraction of total mass (f_mtot).
Fraction of total mass in a stellar generation going into each mass bin (f_m,
f_m_low, f_m_high).  '''
dbin_high = 1.
bins_high = np.arange(m_cutoff, m_max, dbin_high)
n_bins_high = len(bins_high)
m_ave_high = (Gamma / alpha2) * \
             ((bins_high + dbin_high)**alpha2 - bins_high**alpha2) / \
             ((bins_high + dbin_high)**Gamma - bins_high**Gamma)

#####


##### Interpolated Yields #####

##### Minor grid points (interpolate across mass to generate yields at each
##### mass bin of m_ave_low at the original metallicity values)

# yields = 0 for M < 12 Msun and M > 30 Msun

c06_interp_metal = np.zeros((n_bins_high, len(c06_orig)))
for i in xrange(len(c06_orig)):
    itmp = interpolate.InterpolatedUnivariateSpline(m_c06, c06_orig[i], k=1)
    m_tmp = itmp(m_ave_high)
    m_tmp[np.where(m_ave_high < 12.)] = 0.
    m_tmp[np.where(m_ave_high > 30.)] = 0.
    m_tmp[np.where(m_tmp < 0.)] = 0.
    c06_interp_metal[:, i] = m_tmp


# project into 3D arrays (metallicities, masses, isotopes)
z_final = pickle_read(stem_yldgen + 'interp_metallicity.pck')
n_metal_bin = len(z_final)
c06_final = np.ones((n_metal_bin, n_bins_high, len(c06_orig)))
c06_final[:] = c06_interp_metal

##############################

# pickle the interpolated yields array and the metallicity grid used
pickle_write(c06_final, stem_c06 + 'cescutti06_yields.pck')
