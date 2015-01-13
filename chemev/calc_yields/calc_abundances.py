import numpy as np
import copy
import re
import traceback

class Abundances:
    def __init__(self, stem_parent, sym_iso, mgas_iso, weight, timesteps=None,
                 sim_params=None, sim_id=None):
        self.stem_parent = stem_parent
        self.stem_yields = stem_parent + 'yields/'
        self.stem_yldgen = self.stem_yields + 'general/'
        self.isotope = sym_iso
        self.setup()
        self.split_element_mass()
        self.mgas_iso = mgas_iso
        self.n_steps = len(self.mgas_iso)
        self.survivors = weight
        self.t = timesteps
        self.param = sim_params
        self.sim_id = 'ab' + sim_id.lstrip('box')
#        self.apogee_elements()

    def setup(self):
        self.all_atomic_num = np.loadtxt(self.stem_yldgen + 'sym_atomicnum.txt',
                                         dtype=int, usecols=(0,))
        self.all_elements = np.loadtxt(self.stem_yldgen + 'sym_atomicnum.txt',
                                       dtype=str, usecols=(1,))


    def split_element_mass(self):
        '''Takes an array of isotopes (element & mass) and creates a separate
        arrays of element symbols and isotope masses with the same length as
        the isotope array. Also creates a dictionary with the indices of each
        element in the isotope array.'''        
        self.n_isotope = len(self.isotope)
        self.sym = np.array(['' for i in xrange(self.n_isotope)], dtype='|S2')
        self.isotope_mass = np.zeros(self.n_isotope, dtype=int)
        self.elements = []
        for i in xrange(self.n_isotope):
            match = re.match(r"([a-z]+)([0-9]+)", self.isotope[i], re.I)
            if match:
                self.sym[i], self.isotope_mass[i] = match.groups()
            if self.sym[i] not in self.elements:
                self.elements.append(self.sym[i])
        self.elements = np.array(self.elements)
        self.n_elements = len(self.elements)
        self.ind_element = {}
        for item in self.elements:
            self.ind_element[item] = np.where(self.sym == item)[0]


    def load_solar_abund(self, source='lodders'):
        if source == 'lodders':
            fin = self.stem_yldgen + 'lodders03_solar_photosphere.txt'
            self.solar_element = np.loadtxt(fin, dtype=str, usecols=(0,))
            self.solar_abund = np.loadtxt(fin, usecols=(1,))
            self.solar_h = np.zeros(self.n_elements)
            self.solar_fe = np.zeros(self.n_elements)
            for i in xrange(self.n_elements):
                ind = np.where(self.solar_element == self.elements[i])[0]
                ind_fe = np.where(self.elements == 'Fe')
                self.solar_h[i] = self.solar_abund[ind]
                self.solar_fe[i] = np.log10(
                    10.**(self.solar_abund[ind] - 12.) /
                    10.**(self.solar_abund[ind_fe] - 12.))
        #elif source == 'asplund':
        #elif source == 'aspcap':
        #else:
        #    Raise exception
            


    def calc_abundances(self):
        '''Calculate abundances relative to hydrogen and iron.'''        
        self.ngas_iso = np.divide(self.mgas_iso, self.isotope_mass)
        self.niso_h = np.array([
            self.ngas_iso[j] / self.ngas_iso[j, self.ind_element['H']].sum()
            for j in xrange(1, self.n_steps)])
        self.niso_fe = np.array([
            self.ngas_iso[j] / self.ngas_iso[j, self.ind_element['Fe']].sum()
            for j in xrange(1, self.n_steps)])
        self.xh_abs = np.log10([
            np.sum(self.niso_h[:, self.ind_element[item]], axis=1)
            for item in self.elements]) + 12.
        self.xfe_abs = np.log10([
            np.sum(self.niso_fe[:, self.ind_element[item]], axis=1)
            for item in self.elements])
        self.xh_all = np.subtract(self.xh_abs.T, self.solar_h).T
        self.feh = self.xh_all[np.where(self.elements == 'Fe')][0]
        self.xfe_all = np.subtract(self.xfe_abs.T, self.solar_fe).T


    def select_elements(self, el):
        ind = []
        for item in el:
            ind.append(np.where(self.elements == item)[0][0])
        self.xfe = self.xfe_all[ind]
        self.elements_out = self.elements[ind]
        ind2 = []
        for item in el:
            ind2.append(np.where(self.all_elements == item)[0][0])
        self.atomic_num_out = self.all_atomic_num[ind2]
        

    def apogee_solar_abund(self):
        '''Solar Abundances from ASPCAP'''
        self.apogee_el = np.loadtxt(self.stem_yldgen + 'aspcap_solar.txt',
                                    dtype='S2', comments='#', usecols=(0,))
        solar_ab = np.loadtxt(self.stem_yldgen + 'aspcap_solar.txt',
                              comments='#', usecols=(1,))
        self.apogee_solar_H = np.log10(10.**(solar_ab - 12.))
        self.apogee_solar_H = np.delete(self.apogee_solar_H, np.s_[14])
        self.apogee_solar_FeH = np.log10(10.**(apogee_solar_ab[14] - 12.))
        self.apogee_solar_Fe = np.log10(10.**(apogee_solar_ab - 12.) / 
                                        10.**(apogee_solar_ab[14] - 12.))
        self.apogee_solar_Fe = np.delete(self.apogee_solar_Fe, np.s_[14])


# rename attributes as self.apogee_*
    def apogee_elements(self):
        '''Select APOGEE Elements'''
        self.apogee_solar_abund()
        self.ind_apogee = [[12, 13, 14],
                           [15, 16, 17, 18],
                           [19, 20, 21, 22, 23],
                           [32, 33, 34, 35],
                           [36, 37, 38, 39, 40],
                           [41, 42, 43, 44, 45, 46],
                           [47, 48, 49, 50, 51, 52],
                           [59, 60, 61, 62, 63, 64, 65],
                           [78, 79, 80, 81, 82, 83],
                           [84, 85, 86, 87, 88, 89, 90, 91, 92, 93],
                           [103, 104, 105, 106, 107, 108, 109, 110],
                           [111, 112, 113, 114, 115, 116, 117, 118],
                           [119, 120, 121, 122, 123, 124, 125, 126],
                           [127, 128, 129, 130, 131, 132, 133, 134],
                           [145, 146, 147, 148, 149, 150, 151, 152],
                           [153, 154, 155, 156, 157, 158, 159, 160, 161, 162]]
        self.ind_fe = [135, 136, 137, 138, 139, 140, 141, 142, 143, 144]
        self.sym_apogee = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                                    'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co',
                                    'Ni'])
        self.abund_apogee = ['[' + self.sym_apogee[i] + '/Fe]'
                             for i in xrange(len(self.sym_apogee))]
        self.Z_apogee = np.array([6, 7, 8, 11, 12, 13, 14, 16, 19, 20, 22, 23,
                                  24, 25, 27, 28])
        self.n_apogee = len(self.sym_apogee)
