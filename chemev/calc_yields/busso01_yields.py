"""
FILE
    busso01_yields.py

DESCRIPTION

    Generates a finely spaced grid of Ba yields in AGB stars from Busso et
    al. (2001).

    Busso et al. (2001): M = 1.5--3.0 Msun; Z = 2e-4--0.04
    extrapolate the yields down to 1.0 Msun


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
stem_b01 = stem_yields + 'busso01/'
####################


##### Load Data ####
data_in = np.loadtxt(stem_b01 + 'busso01_sprocess.txt', unpack=True)
z_b01 = data_in[0]
xba_b01 = data_in[1:]
m_b01 = np.array([1.5, 3.])
n_b01 = len(z_b01)
b01_ba = xba_b01.T * m_b01 # convert X_Ba to a mass of Ba (and transpose)
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

dbin_low = 0.1
bins_low = np.arange(m_min, m_cutoff, dbin_low)
n_bins_low = len(bins_low)
m_ave_low = (Gamma / alpha2) * \
            ((bins_low + dbin_low)**alpha2 - bins_low**alpha2) / \
            ((bins_low + dbin_low)**Gamma - bins_low**Gamma)
#####


##### Interpolated Yields #####

##### Minor grid points (interpolate across mass to generate yields at each
##### mass bin of m_ave_low at the original metallicity values)

# Linearly extrapolate the yields down to 1.0 Msun

# yields = 0 for M < 1.0 Msun and M > 3.0 Msun

b01_interp_mass = np.zeros((n_b01, n_bins_low))
for i in xrange(n_b01):
    itmp = interpolate.InterpolatedUnivariateSpline(m_b01, b01_ba[i], k=1)
    m_tmp = itmp(m_ave_low)
    m_tmp[np.where(m_ave_low < 1.)] = 0.
    m_tmp[np.where(m_ave_low > 3.)] = 0.
    m_tmp[np.where(m_tmp < 0.)] = 0.
    b01_interp_mass[i] = m_tmp


##### Interpolate across metallicity to generate yields at each mass bin of
##### m_ave_low for N = n_metal_bin metallicity values

# Interpolate Busso et al. (2001) yields onto Limongi & Chieffi (2006)
# metallicity grid, which is evenly sampled in log(metallicity) between each
# metallicity grid point ( 1e-6, 1e-4, 1e-3, 6e-3, 2e-2) for a total of 1001
# values
z_final = pickle_read(stem_yldgen + 'interp_metallicity.pck')
n_metal_bin = len(z_final)

# Extend the metallicity grid up to Z = 0.04 (twice solar), since Busso et
# al. (2001) provide yields up to this metallicity
n_metal_bin2 = 401
z_grid2 = np.array([2e-2, 3e-2, 4e-2])
logz_grid2 = np.log10(z_grid2)

# evenly sample metallicity (in log Z) between grid points
logz_final2 = np.zeros(n_metal_bin2)
dind = (n_metal_bin2 - 1) / (len(z_grid2) - 1)
for i in xrange(len(z_grid2) - 1):
    dlogz2 = (logz_grid2[i+1] - logz_grid2[i]) / \
             ((n_metal_bin2 - 1) / (len(z_grid2) - 1))
    logz_final2[i*dind:i*dind+dind+1] = np.arange(logz_grid2[i],
                                                  logz_grid2[i+1]+1e-9, dlogz2)

# metallicity of extended grid
z_final2 = 10.**logz_final2
# metallicity of original plus extended grid
z_final3 = np.concatenate((z_final, z_final2[1:]))
n_metal_bin3 = len(z_final3)

# interpolated Busso et al. (2001) yields
# metallicity grid of Limongi & Chieffi (2006) (up to solar)
b01_final = np.zeros((n_metal_bin, n_bins_low))
# metallicity grid of Busso et al. (2001) (up to twice solar)
b01_final_ext = np.zeros((n_metal_bin3, n_bins_low))
# at each mass, interpolate each element for each metallicity
for i in xrange(n_bins_low):
    itmp = interpolate.InterpolatedUnivariateSpline(
        z_b01, b01_interp_mass[:, i], k=1)
    b01_final[:, i] = itmp(z_final)
    b01_final[np.where(z_final < z_b01[0])] = 0.
    b01_final[np.where(b01_final < 0.)] = 0.
    itmp2 = interpolate.InterpolatedUnivariateSpline(
        z_b01, b01_interp_mass[:, i], k=1)
    b01_final_ext[:, i] = itmp2(z_final3)
    b01_final_ext[np.where(z_final3 < z_b01[0])] = 0.
    b01_final_ext[np.where(b01_final_ext < 0.)] = 0.


# project into 3D arrays (metallicities, masses, isotopes)
b01_final = np.array([b01_final]).transpose(1, 2, 0)
b01_final_ext = np.array([b01_final_ext]).transpose(1, 2, 0)
##############################

# pickle the interpolated yields array and the metallicity grid used
pickle_write(b01_final, stem_b01 + 'busso01_yields.pck')
pickle_write(b01_final_ext, stem_b01 + 'busso01ext_yields.pck')
pickle_write(z_final3, stem_b01 + 'busso01ext_metallicity.pck')
