"""
FILE
    ww95_yields.py

DESCRIPTION

    Generates a finely spaced grid of SN II isotopic yields from Woosley &
    Weaver (1995), AGB isotopic yields from Renzini & Voli (1981), and SNIa
    yields from Thielemann, Nomoto, & Yokoi (1986).
    
    Woosley & Weaver (1995): M = 11--40 Msun; Z = 0--solar
    Renzini & Voli (1981): M = 1--8 Msun; Z = 0--solar
    Thielemann et al. (1986): W7 model from Nomoto et al. (1984)

    Timmes already converted Ni56 to Fe56 in the maltov1.orig file (WW95
    doesn't account for its decay).

"""

import numpy as np
import os
import pickle
import copy
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
stem_ww95 = stem_yields + 'ww95/'
stem_ww95_orig = stem_ww95 + 'orig/'
stem_ww95_half_fe = stem_ww95 + 'half_fe/'
#stem_ww95_half_fe_only = stem_ww95 + 'half_fe_only/'
#stem_rv81 = stem_flexce + 'yields/renzini81/'
#stem_tny86 = stem_flexce + 'yields/thielemann86/'
####################

##### WW95 Yields #####
z = open(stem_ww95 + 'maltov1.orig', 'r')
sym = [] # symbol names
sym_metallicity = [] # [symbol--metallicity pairs]
bbmf = [] # big bang mass fraction
sneIa_orig = []
ww95_orig = [] # ww95_orig[80 symbols + 5 metallicities/sym = 400][[25 masses], [25 yields]]
tmp_Ia = 0; tmp = 0
for row in z:
    if 'symbol name' in row:
        sym_tmp = row.split()[0]
        sym.append(sym_tmp)
    if 'big bang mass fraction' in row:
        bbmf.append(float(row.split()[0]))
    if 'w7 tny86' in row:
        yields_Ia = []
        tmp_Ia = 6
    if tmp_Ia > 0:
        yields_Ia.append(float(row.split()[0]))
        tmp_Ia -= 1
        if tmp_Ia ==0:
            sneIa_orig.append(np.array(yields_Ia))
    if '* metallicity' in row:
        metal_tmp = float(row.split()[0])
        sym_metallicity.append([sym_tmp, metal_tmp])
    if 'rv81 stellar mass & yield' in row:
        mass = []
        yields = []
        tmp = 25
    if tmp > 0:
        mass.append(float(row.split()[0]))
        yields.append(float(row.split()[1]))
        tmp -= 1
        if tmp == 0:
            ww95_orig.append([np.array(mass), np.array(yields)])

z.close()
sym = np.array(sym)
sym_mass = np.array([int(sym[i][-1]) if i < 7 else int(sym[i][-2:]) \
                     for i in xrange(len(sym) - 1)])
sym_metallicity = np.array(sym_metallicity)
bbmf = np.array(bbmf)[:-1]
sneIa_orig = np.array(sneIa_orig)
tnyIa = sneIa_orig[:, 0]
ww95_orig = np.array(ww95_orig)
# all symbols have 25 masses and yields and 5 metallicity values:
ww95_mass = ww95_orig[0][0]
ww95_mass2 = np.concatenate([ww95_mass for i in xrange(5)])
ww95_metal = np.array([0.00e+00, 1.90e-06, 1.90e-04, 1.90e-03, 1.90e-02])
ww95_metal2 = np.concatenate([np.ones(25) * ww95_metal[i] for i in xrange(5)])
n_sym = len(sym)
n_iso = len(sym) - 1
n_metal = len(sym) - 5
n_yield = len(ww95_orig)
#######################


##### CL04 Data ####
species = np.loadtxt(stem_yldgen + 'species.txt', dtype=str, usecols=(1,))
n_species = len(species)

# match isotopes from WW95 yields to CL04 yields
sym2 = np.array([sym[i].title() for i in xrange(len(sym))])
ind_sp = []
for i in xrange(n_sym):
    if sym2[i] in species:
        tmp = np.where(sym2[i] == species)[0][0]
        ind_sp.append(tmp)
    else:
        pass
        #print 'sym[%i]' % (i), '(%s)' % (sym[i]), 'not in species array'

ind_sp = np.array(ind_sp)

# solar abundance of metals---needed to subtract the initial metal abundances
# of the stellar models (also assume Y = 0.285)---in relative amounts (not
# Msun), that is, sum(solar_ab) = 1.
solar_iso = np.loadtxt(stem_yldgen + 'Solar_isotopes.txt', skiprows=1, dtype=str,
                       usecols=(0,))
solar_ab = np.loadtxt(stem_yldgen + 'Solar_isotopes.txt', skiprows=1,
                      usecols=(1,))

# indices within "species" array of the elements for which CL04 give a solar
# abundance (Note: WW95 also used the Anders & Grevesse 1989 solar abundance)
ind_iso = []
for i in xrange(len(solar_iso)):
    ind_iso.append(np.where(solar_iso[i] == species)[0][0])

ind_iso = np.array(ind_iso)
####################


### Calculate Net Yields ###

# WW95 absolute yields (125 mass/metallicity pairs, 293 isotopes)
ww95_orig2 = ww95_orig.reshape(80, 5, 2, 25)
ww95_orig3 = ww95_orig2[:, :, 1]
ww95_orig4 = ww95_orig3.reshape(80, 125).T
ww95_abs = np.zeros((125, n_species))
for i in xrange(125):
    for j in xrange(79):
        ww95_abs[i, ind_sp[j]] = ww95_orig4[i, j]

# WW95 mass ejected
ww95_mej = np.sum(ww95_abs, axis=1)

# WW95 remnant mass
ww95_rem = ww95_mass2 - ww95_mej

# The remnant masses reported by WW95 but the sum(abs yields) + remnant mass !=
# mass of star, so for accouting purposes it will be best to calculate remnant
# mass = mass of star - sum(abs yields).
# WW95 reported remnant masses:
# ww95_rem = ww95_orig4[:, -1]

# WW95 initial composition
ww95_init_comp = np.zeros(ww95_abs.shape)

for i in xrange(5):
    indt = np.arange(25*i, 25*i+25)
    if i == 0:
        ww95_init_comp[indt, 0] = (1. - 0.23) * ww95_mass # H
        ww95_init_comp[indt, 4] = 0.23 * ww95_mass # He
    else:
        ztmp = ww95_metal2[indt][0]
        ww95_init_comp[indt, 0] = (1. - 0.285 - ztmp) * ww95_mass # H
        ww95_init_comp[indt, 4] = 0.285 * ww95_mass # He
        for j in xrange(len(ind_iso)):
            ww95_init_comp[indt, ind_iso[j]] = (ztmp * solar_ab[j] *
                                                ww95_mass) #C-->Mo

# WW95 net yields = absolute yields - initial composition of stellar model
ww95_net = ww95_abs - ww95_init_comp

# WW95 SNII (11--40 Msun) net yields, mass ejected, and remnant mass for Z > 0
ww95_sn_net = np.zeros((48, n_species))
ww95_sn_mej = np.zeros(48)
ww95_sn_rem = np.zeros(48)
for i in xrange(1, 5):
    ind1 = np.arange(25*i+13, 25*i+25)
    ind2 = np.arange(12*(i-1), 12*(i-1)+12)
    ww95_sn_net[ind2] = ww95_net[ind1]
    ww95_sn_mej[ind2] = ww95_mej[ind1]
    ww95_sn_rem[ind2] = ww95_rem[ind1]


# Renzini & Voli (1981) AGB net yields (1--8 Msun) for Z > 0
rv81_net = np.zeros((44, n_species))
rv81_mej = np.zeros(44)
rv81_rem = np.zeros(44)
for i in xrange(1, 5):
    ind1 = np.arange(25*i, 25*i+11)
    ind2 = np.arange(11*(i-1), 11*(i-1)+11)
    rv81_net[ind2] = ww95_net[ind1]
    rv81_mej[ind2] = ww95_mej[ind1]
    rv81_rem[ind2] = ww95_rem[ind1]

############################


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

dbin_high = 1.
bins_high = np.arange(m_cutoff, m_max, dbin_high)
n_bins_high = len(bins_high)
m_ave_high = (Gamma / alpha2) * \
             ((bins_high + dbin_high)**alpha2 - bins_high**alpha2) / \
             ((bins_high + dbin_high)**Gamma - bins_high**Gamma)

m_ave = np.append(m_ave_low, m_ave_high)
n_bins = n_bins_low + n_bins_high
#####


##### Interpolated Yields #####

##### Minor grid points (mass bins spaced in ~1 Msun, but at the original 4
##### non-zero metallicity values [Z = 1.9e-6, 1.9e-4, 1.9e-3, 1.9e-2])

# Interpolate across mass to generate yields at each mass bin of m_ave_high for
# 4 metallicity values: Z = 1.9e-6, 1.9e-4, 1.9e-3, & 1.9e-2 (almost solar)

# Linearly extrapolate the mass ejected and the remnant mass at fixed
# metallicity

# Do NOT extrapolate (linearly or otherwise) net yields.

ww95_sn_m = np.array([11.065, 12.065, 13.071, 15.081, 18.098, 19., 20.109,
                      22.119, 25.136, 30.163, 35.19, 40.217])
rv81_m = np.array([1., 1.5, 1.75, 2., 2.5, 3., 4., 5., 6., 7., 8.])

ind_ww95_sn_net = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                   [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                   [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
                   [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]]

ind_rv81_net = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
                [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
                [33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]]

# WW95 yields
ww95_interp_mass = np.zeros((4, n_bins_high, n_species))
for i in xrange(4):
    m_tmp = np.zeros((n_bins_high, n_species))
    for k in xrange(n_species):
        itmp = interpolate.InterpolatedUnivariateSpline(
            ww95_sn_m, ww95_sn_net[ind_ww95_sn_net[i], k], k=1)
        m_tmp[:, k] = itmp(m_ave_high)
        m_tmp[np.where(m_ave_high < ww95_sn_m[0]), k] = itmp(ww95_sn_m[0])
        m_tmp[np.where(m_ave_high > ww95_sn_m[-1]), k] = itmp(ww95_sn_m[-1])
    ww95_interp_mass[i] = m_tmp

# WW95 mass ejected
ww95_interp_mej = np.zeros((4, n_bins_high))
for i in xrange(4):
    itmp = interpolate.InterpolatedUnivariateSpline(
        ww95_sn_m, ww95_sn_mej[ind_ww95_sn_net[i]], k=1)
    ww95_interp_mej[i] = itmp(m_ave_high)

# WW95 remnant mass
ww95_interp_rem = np.zeros((4, n_bins_high))
for i in xrange(4):
    ww95_interp_rem[i] = m_ave_high - ww95_interp_mej[i]

# rv81 yields
rv81_interp_mass = np.zeros((4, n_bins_low, n_species))
for i in xrange(4):
    m_tmp = np.zeros((n_bins_low, n_species))
    for k in xrange(n_species):
        itmp = interpolate.InterpolatedUnivariateSpline(
            rv81_m, rv81_net[ind_rv81_net[i], k], k=1)
        m_tmp[:, k] = itmp(m_ave_low)
        m_tmp[np.where(m_ave_low < rv81_m[0]), k] = itmp(rv81_m[0])
    rv81_interp_mass[i] = m_tmp

# RV81 mass ejected
rv81_interp_mej = np.zeros((4, n_bins_low))
for i in xrange(4):
    itmp = interpolate.InterpolatedUnivariateSpline(
        rv81_m, rv81_mej[ind_rv81_net[i]], k=1)
    rv81_interp_mej[i] = itmp(m_ave_low)
    rv81_interp_mej[i][np.where(rv81_interp_mej[i] < 0.)] = 0.

# RV81 remnant mass
rv81_interp_rem = np.zeros((4, n_bins_low))
for i in xrange(4):
    rv81_interp_rem[i] = m_ave_low - rv81_interp_mej[i]


##### Interpolate across metallicity to generate yields at each mass bin of
##### m_ave_high for N = n_metal_bin metallicity values

# Interpolate WW95 yields onto Limongi & Chieffi (2006) metallicity grid, which
# is evenly sampled in log(metallicity) between each metallicity grid point (
# 1e-6, 1e-4, 1e-3, 6e-3, 2e-2) for a total of 1001 values
z_final = pickle_read(stem_yldgen + 'interp_metallicity.pck')
n_metal_bin = len(z_final)


# interpolated WW95 SNII yields
ww95_final = np.zeros((n_metal_bin, n_bins_high, n_species))
# at each mass, interpolate each element for each metallicity
for i in xrange(n_bins_high):
    for j in xrange(n_species):
        itmp = interpolate.InterpolatedUnivariateSpline(
            ww95_metal[1:], ww95_interp_mass[:, i, j], k=1)
        ww95_final[:, i, j] = itmp(z_final)

# interpolated WW95 SNII mass ejected
ww95_final_mej = np.zeros((n_metal_bin, n_bins_high))
for i in xrange(n_bins_high):
    itmp = interpolate.InterpolatedUnivariateSpline(ww95_metal[1:],
                                                    ww95_interp_mej[:, i], k=1)
    ww95_final_mej[:, i] = itmp(z_final)

# interpolated WW95 SNII remnant mass
ww95_final_rem = np.zeros((n_metal_bin, n_bins_high))
for i in xrange(n_metal_bin):
    ww95_final_rem[i] = m_ave_high - ww95_final_mej[i]


# interpolated RV81 AGB yields
rv81_final = np.zeros((n_metal_bin, n_bins_low, n_species))
# at each mass, interpolate each element for each metallicity
for i in xrange(n_bins_low):
    for j in xrange(n_species):
        itmp = interpolate.InterpolatedUnivariateSpline(
            ww95_metal[1:], rv81_interp_mass[:, i, j], k=1)
        rv81_final[:, i, j] = itmp(z_final)

# interpolated RV81 AGB mass ejected
rv81_final_mej = np.zeros((n_metal_bin, n_bins_low))
for i in xrange(n_bins_low):
    itmp = interpolate.InterpolatedUnivariateSpline(ww95_metal[1:],
                                                    rv81_interp_mej[:, i], k=1)
    rv81_final_mej[:, i] = itmp(z_final)

# interpolated RV81 SNII remnant mass
rv81_final_rem = np.zeros((n_metal_bin, n_bins_low))
for i in xrange(n_metal_bin):
    rv81_final_rem[i] = m_ave_low - rv81_final_mej[i]

#####################################


### WW95 yields with 1/2 Fe and Fe-peak element (Cr, Mn, Fe, Co, Ni, Cu, Zn) yields
ww95_final_half_fe = copy.deepcopy(ww95_final)
ww95_final_half_fe[:, :, 119:185] = ww95_final_half_fe[:, :, 119:185] / 2.
ww95_final_mej_half_fe = (ww95_final_mej -
                          np.sum(ww95_final_half_fe[:, :, 119:185], axis=2))
ww95_final_rem_half_fe = (ww95_final_rem +
                          np.sum(ww95_final_half_fe[:, :, 119:185], axis=2))

### WW95 yields with 1/2 Fe yields
ww95_final_half_fe_only = copy.deepcopy(ww95_final)
ww95_final_half_fe_only[:, :, 135:145] = ww95_final_half_fe_only[:, :, 135:145] / 2.
ww95_final_mej_half_fe_only = (ww95_final_mej -
                               np.sum(ww95_final_half_fe_only[:, :, 135:145], axis=2))
ww95_final_rem_half_fe_only = (ww95_final_rem +
                               np.sum(ww95_final_half_fe_only[:, :, 135:145], axis=2))


### Thielemann, Nomoto, & Yokoi (1986) detailed nucleosynthesis calculations of
### the Nomoto et al. (1984) W7 SN Ia model
tny86_final = np.zeros(n_species)
tny86_final[ind_sp] = tnyIa

#################################




# pickle the interpolated yields array and the metallicity grid used
pickle_write(ww95_final, stem_ww95_orig + 'interp_yields.pck')
pickle_write(ww95_final_mej, stem_ww95_orig + 'interp_meject.pck')
pickle_write(ww95_final_rem, stem_ww95_orig + 'interp_mremnant.pck')
#pickle_write(rv81_final, stem_rv81 + 'interp_yields.pck')
#pickle_write(rv81_final_mej, stem_rv81 + 'interp_meject.pck')
#pickle_write(rv81_final_rem, stem_rv81 + 'interp_mremnant.pck')
#pickle_write(tny86_final, stem_tny86 + 'w7_yields.pck')

pickle_write(ww95_final_half_fe, stem_ww95_half_fe + 'interp_yields.pck')
pickle_write(ww95_final_mej_half_fe, stem_ww95_half_fe + 'interp_meject.pck')
pickle_write(ww95_final_rem_half_fe, stem_ww95_half_fe + 'interp_mremnant.pck')

#pickle_write(ww95_final_half_fe_only, stem_ww95_half_fe_only + 'interp_yields.pck')
#pickle_write(ww95_final_mej_half_fe_only, stem_ww95_half_fe_only + 'interp_meject.pck')
#pickle_write(ww95_final_rem_half_fe_only, stem_ww95_half_fe_only + 'interp_mremnant.pck')




# ##### Interpolate WW95 Yields Between Mass and Metallicity Values #####
# 
# '''Interpolate the WW95 yields between the bins of stellar mass for massive
# stars (92 bins).  ww95_interp_mass[80 symbols + 5 metallicities/sym = 400][[92
# mass--yield pairs]].  '''
# 
# ww95_tmp = []
# for i in xrange(n_yield):
#     interp_mass = interpolate.interp1d(ww95_orig[i][0], ww95_orig[i][1] , \
#                                        bounds_error=False, fill_value=0.)
#     ww95_tmp.append([m_ave, interp_mass(m_ave)])
# 
# ww95_interp_mass = np.array(ww95_tmp).transpose((0, 2, 1))
# ww95a = ww95_interp_mass.reshape(80, 5, n_bins, 2)
# ww95b = ww95a.transpose(0, 2, 3, 1)
# 
# def interp_ww95(metallicity):    
#     """Interpolate the WW95 yields between the metallicity values that they
#     give.  If metallicity falls outside of the metallicity value grid, then
#     return the yield for the maximum metallicity.
#     """
#     ww95_interp_metallicity = []# ww95_interp_metallicity[symbol][interpolated yield]
#     for i in xrange(n_sym):
#         stmp = []
#         for j in xrange(len(m_ave)):
#             interp_metallicity = interpolate.InterpolatedUnivariateSpline( \
#                 sn2_metallicity, ww95a[i, :, j, 1], k=1)
#             stmp.append(interp_metallicity(metallicity))
#         ww95_interp_metallicity.append(stmp)
#     ww95_interp_metallicity = np.array(ww95_interp_metallicity)
#     return ww95_interp_metallicity
# 
# #########################################################################
