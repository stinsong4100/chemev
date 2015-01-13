import numpy as np
import copy
import pickle

def pickle_read(filename):
    '''read pickled object'''
    fin = open(filename, 'r')
    obj = pickle.load(fin)
    fin.close()
    return obj
    
class Yields:
    def __init__(self, stem_parent, snii_dir, agb_dir, snia_dir,
                 rprocess_dir=None, sprocess_dir=None,
                 mbins=np.concatenate((np.arange(0.1, 8., 0.1),
                                       np.arange(8., 100.1, 1.)))):
        self.stem_parent = stem_parent
        self.stem_yields = stem_parent + 'yields/'
        self.stem_yldgen = self.stem_yields + 'general/'
        self.stem_snii = self.stem_yields + snii_dir
        self.stem_agb = self.stem_yields + agb_dir
        self.stem_snia = self.stem_yields + snia_dir
        self.sources = dict(snii=snii_dir.split('/')[0],
                            snia=snia_dir.split('/')[0],
                            agb=agb_dir.split('/')[0])
        if rprocess_dir is not None:
            self.stem_rprocess = self.stem_yields + rprocess_dir
            self.sources['rprocess'] = rprocess_dir.split('/')[0]
        if sprocess_dir is not None:
            self.stem_sprocess = self.stem_yields + sprocess_dir
            self.sources['sprocess'] = sprocess_dir.split('/')[0]
        self.mass_bins = mbins
        self.n_bins = len(mbins) - 1
        self.n_bins_high = len(np.where(self.mass_bins >= 8.)[0]) - 1
        self.n_bins_low = len(np.where(self.mass_bins < 8.)[0])
        self.ind8 = np.where(self.mass_bins == 8.)[0][0]
        self.mlow = mbins[0]
        self.load_sym()
        self.load_bbmf()


    def load_sym(self):
        '''Load isotopic and elemental symbols and masses'''
        self.atomic_num = np.loadtxt(self.stem_yldgen + 'sym_atomicnum.txt',
                                     dtype=int, usecols=(0,))
        self.element_all = np.loadtxt(self.stem_yldgen + 'sym_atomicnum.txt',
                                  dtype=str, usecols=(1,))
        self.snii_sym = np.loadtxt(self.stem_yldgen + 'species.txt', dtype=str,
                                 usecols=(1,))
        self.snii_sym_mass = np.loadtxt(self.stem_yldgen + 'species.txt',
                                        dtype=int, usecols=(2,))
        self.n_snii_sym = len(self.snii_sym)
        u, indices = np.unique([item.rstrip('0123456789')
                                for item in self.snii_sym], return_index=True)
        indices_s = np.argsort(indices)
        self.element = np.delete(u[indices_s], np.s_[13, 14])
        self.n_elements = len(self.element)


    def load_bbmf(self):
        ''' big bang mass fraction (originally from WW95/Timmies data file)'''
        self.bbmf = pickle_read(self.stem_yldgen + 'bbmf.pck')


    def load_solar_abund(self, source='lodders'):
        if source == 'lodders':
            fin = self.stem_yldgen + 'lodders03_solar_photosphere.txt'
            self.solar_element = np.loadtxt(fin, dtype=str, usecols=(0,))
            self.solar_abund = np.loadtxt(fin, usecols=(1,))
            self.solar_h = np.zeros(self.n_elements)
            self.solar_fe = np.zeros(self.n_elements)
            for i in xrange(self.n_elements):
                ind = np.where(self.solar_element == self.element[i])[0]
                ind_fe = np.where(self.element == 'Fe')
                self.solar_h[i] = self.solar_abund[ind]
                self.solar_fe[i] = np.log10(
                    10.**(self.solar_abund[ind] - 12.) /
                    10.**(self.solar_abund[ind_fe] - 12.))
        self.calc_solar_mass_fraction()
        #elif source == 'asplund':
        #elif source == 'aspcap':
        #else:
        #    Raise exception


    def calc_solar_mass_fraction(self):
        dominant_isotope = np.array([
            'H1', 'He4', 'Li7', 'Be9', 'B11', 'C12', 'N14', 'O16', 'F19',
            'Ne20', 'Na23', 'Mg24', 'Al27', 'Si28', 'P31', 'S32', 'Cl35',
            'Ar36', 'K39', 'Ca40', 'Sc45', 'Ti48', 'V51', 'Cr52', 'Mn55',
            'Fe56', 'Co59', 'Ni58', 'Cu63', 'Zn64', 'Ga69', 'Ge74', 'As75',
            'Se80', 'Br81', 'Kr84', 'Rb85', 'Sr88', 'Y89', 'Zr90', 'Nb93',
            'Mo98', 'Ba138', 'Eu153'], dtype='|S5')
        miso = np.array([item.lstrip(''.join(['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
                                              'abcdefghijklmnopqrstuvwxyz']))
                         for item in dominant_isotope])
        miso = miso.astype(int)
        nh = 10.**(self.solar_h - 12.)
        mh = nh * miso
        mfrac = mh / mh.sum()
        self.solar_mfrac = np.zeros(self.n_sym)
        for i in xrange(len(dominant_isotope)):
            ind = np.where(self.sym == dominant_isotope[i])[0]
            self.solar_mfrac[ind] = mfrac[i]


    def load_snii_yields(self):
        ''' Limongi & Chieffi SN II Yields: 1001 metallicities, 92 masses from
        8 to 100 msun, 293 isotopes.'''
        self.snii_z = pickle_read(self.stem_yldgen + 'interp_metallicity.pck')
        self.snii_yields = pickle_read(self.stem_snii + 'interp_yields.pck')
        self.snii_mej = pickle_read(self.stem_snii + 'interp_meject.pck')
        self.snii_rem = pickle_read(self.stem_snii + 'interp_mremnant.pck')
        self.snii_agb_z = copy.deepcopy(self.snii_z)
        self.n_z = len(self.snii_agb_z)
        # adjust for a different upper mass limit to the IMF
        self.snii_yields = self.snii_yields[:, :self.n_bins_high, :]
        self.snii_mej = self.snii_mej[:, :self.n_bins_high]
        self.snii_rem = self.snii_rem[:, :self.n_bins_high]



    def load_agb_yields(self):
        '''Karakas AGB Yields: 1001 metallicities, 79 masses from 0.1 to 8
        msun, 70 isotopes.

        Only works for Karakas et al. (2010) AGB yields.

        B8, O14, Ne19, and Si33 are not in the LC yields.  The latter three K10
        yields are always 0 and the B8 yield is essentially zero (<1e-17), so I
        have removed those from agb_sym, agb_yields_in, and agb_yields.  '''
        
        self.agb_sym = np.loadtxt(self.stem_agb + 'species.txt', dtype=str,
                                  usecols=(1,))
        self.agb_z = pickle_read(self.stem_yldgen + 'interp_metallicity.pck')
        agb_yields_in = pickle_read(self.stem_agb + 'interp_yields.pck')
        # remove isotopes not in the LC06 yields
        agb_yields_in = np.delete(agb_yields_in, np.s_[6, 13, 23, 47], axis=2)
        self.agb_mej = pickle_read(self.stem_agb + 'interp_meject.pck')
        self.agb_rem = pickle_read(self.stem_agb + 'interp_mremnant.pck')
        # Create an array with the same elements as the Limongi & Chieffi
        # (2006) SNII yields.
        self.agb_yields = np.zeros((agb_yields_in.shape[0],
                                    agb_yields_in.shape[1],
                                    self.snii_yields.shape[2]))
        ind_agb = np.array([], dtype=int)
        for i in xrange(len(self.agb_sym)):
            tmp = np.where(self.agb_sym[i] == self.snii_sym)[0]
            ind_agb = np.append(ind_agb, tmp)
        self.agb_yields[:, :, ind_agb] = copy.deepcopy(agb_yields_in)
        if self.mlow == 0.2:
            self.agb_mej = self.agb_mej[:, 1:]
            self.agb_rem = self.agb_rem[:, 1:]
            self.agb_yields = self.agb_yields[:, 1:]


    def load_rprocess_yields(self):
        '''Cescutti et al. (2006) r-process Ba & Eu yields for M = 12, 15, 30
        Msun that are metallicity independent.  '''
        self.rprocess_yields = pickle_read(self.stem_rprocess +
                                           'cescutti06_yields.pck')


    def load_sprocess_yields(self, supersolar=False):
        '''Busso et al. (2001) s-process Ba yields for M < 8 Msun & Z =
        2e-4--4e-2 (twice solar and twice the maximum Limongi & Chieffi
        metallicity).

        supersolar: if True, use the yields that extend from solar (Z = 2e-2)
        to twice solar (Z = 4e-2), otherwise only use the yields up to solar.
        '''        
        if supersolar:
            self.sprocess_z = pickle_read(self.stem_sprocess +
                                          'busso01ext_metallicity.pck')
            self.sprocess_yields = pickle_read(self.stem_sprocess +
                                               'busso01ext_yields.pck')
        else:
            self.sprocess_yields = pickle_read(self.stem_sprocess +
                                               'busso01_yields.pck')


    def load_snia_yields(self, model='merger11'):
        '''Pakmor SN Ia models:        
        subch06: sub-Chandrasekhar 0.6 Msun model
        merger11: double-degenerate scenario merger 1.1 Msun model
        deldet07: single-degenerate scenario delayed detonation 0.7 Msun model
        subchnew: sub-Chandrasekhar model (not sure what "new" refers to)
        subchnew_0.1: sub-Chandrasekhar model (not sure what "new" or "0.1"
            refer to)
        merger11_0.1: double-degenerate scenario merger 1.1 Msun model (not
            sure what "new" refers to)

        Iwamoto et al. (1999) models:
        W7: single-degenerate model from Nomoto et al. (1984)
        W70: zero-metallicity version of W7
        '''
        self.snia_yields = pickle_read(self.stem_snia + model + '_yields.pck')


    def concat_ncapture_yields(self, r_elements, s_elements):
        '''Create an array of r- and s-process isotopic yields.
        '''
        nclib = pickle_read(self.stem_yldgen + 'sneden08.pck')
        # unique elements arranged by atomic number
        elements = np.unique(r_elements + s_elements)
        at_num = []
        for item in elements:
            at_num.append(nclib[item]['Z'])
        at_num = np.array(at_num)
        elements = elements[np.argsort(at_num)]
        self.nc_sym = []
        self.nc_sym_mass = []
        for item in elements:
            n_tmp = len(nclib[item]['Isotope'])
            for i in xrange(n_tmp):
                self.nc_sym.append(item + str(nclib[item]['Isotope'][i]))
                self.nc_sym_mass.append(nclib[item]['Isotope'][i])
        self.nc_sym = np.array(self.nc_sym)
        self.nc_sym_mass = np.array(self.nc_sym_mass)
        self.n_nc_sym = len(self.nc_sym)
        u, indices = np.unique([item.rstrip('0123456789')
                                for item in self.nc_sym], return_index=True)
        indices_s = np.argsort(indices)
        self.nc_element = u[indices_s]
        # project elemental yields onto relative isotopic abundances
        self.nc_yields = np.zeros((self.n_z, self.n_bins, self.n_nc_sym))
        cnt = 0
        for i in xrange(len(elements)):
            el = elements[i]
            el_iso = len(nclib[el]['Isotope'])
            if el in r_elements:
                j = np.where(np.array(r_elements) == el)[0]
                self.nc_yields[:, -self.n_bins_high:, cnt:cnt+el_iso] = \
                               (np.ones((self.n_z, self.n_bins_high, el_iso)) *
                               self.rprocess_yields[:, -self.n_bins_high:, j] *
                                nclib[el]['isotopic_fraction[r]'])
            if el in s_elements:
                j = np.where(np.array(s_elements) == el)[0]
                self.nc_yields[:, :self.n_bins_low, cnt:cnt+el_iso] = \
                                (np.ones((self.n_z, self.n_bins_low, el_iso)) *
                                 self.sprocess_yields[:, :, j] *
                                 nclib[el]['isotopic_fraction[s]'])
            cnt += el_iso
        # update arrays
        self.sym = np.append(self.snii_sym, self.nc_sym)
        self.sym_mass = np.append(self.snii_sym_mass, self.nc_sym_mass)
        self.n_sym = len(self.sym)
        self.element = np.append(self.element, self.nc_element)
        self.n_elements = len(self.element)
        self.bbmf = np.append(self.bbmf, np.zeros(self.n_nc_sym))
        self.snia_yields = np.append(self.snia_yields, np.zeros(self.n_nc_sym))
        self.snii_yields = np.append(self.snii_yields,
                                     self.nc_yields[:, self.ind8:], axis=2)
        self.agb_yields = np.append(self.agb_yields,
                                    self.nc_yields[:, :self.ind8], axis=2)
        self.snii_agb_rem = np.concatenate((self.agb_rem, self.snii_rem),
                                           axis=1)
