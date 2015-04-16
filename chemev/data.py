import pickle, numpy as np

# Read in chemical evolution tables
global agb_dat , iagbAgeSort, AGBageIndices, iagbZsort, AGBZindices
try:
    agb_dat = pickle.load(open('agb_enrich.pck'))
except:
    print "Couldn't open agb_enrich.pck"
    pass
iagbAgeSort = np.argsort(agb_dat['times'][0])
AGBageIndices = np.arange(len(agb_dat['times'][0]))
iagbZsort = np.argsort(agb_dat['Zs'])
AGBZindices = np.arange(len(agb_dat['Zs']))

global snii_dat, isnAgeSort, SNageIndices, isnZsort, SNZindices
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
