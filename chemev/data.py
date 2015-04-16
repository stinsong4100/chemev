import pickle, numpy as np

# Read in chemical evolution tables
global agb_dat 
global snii_dat, isnAgeSort, SNageIndices, isnZsort, SNZindices
try:
    agb_dat = pickle.load(open('agb_enrich.pck'))
except:
    print "Couldn't open agb_enrich.pck"
    pass

try:
    snii_dat = pickle.load(open('snii_enrich.pck'))
except:
    print "Couldn't open snii_enrich.pck"
    pass
isnAgeSort = np.argsort(snii_dat['times'][0])
SNageIndices = np.arange(len(snii_dat['times'][0]))
isnZsort = np.argsort(snii_dat['Zs'])
SNZindices = np.arange(len(snii_dat['Zs']))

global star_type
star_type = np.dtype({'names':['tform','init_mass','mass','Z'],
                      'formats':['f','f','f','f']})
