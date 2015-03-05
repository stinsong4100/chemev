"""
FILE
    sneden08.py

DESCRIPTION

    Read in the isotopic abundances of neutron-capture isotopes from Sneden et
    al. (2008) and create a pickled dictionary.
    
"""

import numpy as np
import os
import pickle

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
####################

##### Load Data #####
s08 = {}
fin = open(stem_yldgen + 'sneden08.txt', 'r')
for line in fin:
    if line.split()[0] == '#':
        pass
    else:
        cols = line.strip().split()
        if len(cols) == 5:
            name, at_num, at_mass, ns, nr = cols
            s08[name] = {'element': name, 'Z': int(at_num),
                'Isotope': [int(at_mass)],'N[s]': [float(ns)],
                'N[r]': [float(nr)]}
        else:
            at_mass, ns, nr = cols
            s08[name]['Isotope'].append(int(at_mass))
            s08[name]['N[s]'].append(float(ns))
            s08[name]['N[r]'].append(float(nr))

fin.close()
####################


for e in s08.iterkeys():
    for k in ['N[r]', 'N[s]', 'Isotope']:
        s08[e][k] = np.array(s08[e][k])


# Calculate the fraction of each element produced in the r- and s-processes
# ('fraction[r]' and 'fraction[s]') and the fraction of each isotope that
# comes from the r- or s-processes relative to the total amount of the element
# produced in the r- or s-processes('isotopic_fraction[r]' and
# 'isotopic_fraction[s]')
for e in s08.iterkeys():
    ncap_tot = np.sum(s08[e]['N[r]'].sum() + s08[e]['N[s]'].sum())
    s08[e]['fraction[r]'] = s08[e]['N[r]'].sum() / ncap_tot
    s08[e]['fraction[s]'] = s08[e]['N[s]'].sum() / ncap_tot
    if s08[e]['N[r]'].sum() > 0.:
        s08[e]['isotopic_fraction[r]'] = s08[e]['N[r]'] / s08[e]['N[r]'].sum()
    else:
        s08[e]['isotopic_fraction[r]'] = np.zeros(len(s08[e]['N[r]']))
    if s08[e]['N[s]'].sum() > 0.:
        s08[e]['isotopic_fraction[s]'] = s08[e]['N[s]'] / s08[e]['N[s]'].sum()
    else:
        s08[e]['isotopic_fraction[s]'] = np.zeros(len(s08[e]['N[s]']))


pickle_write(s08, stem_yldgen + 'sneden08.pck')
