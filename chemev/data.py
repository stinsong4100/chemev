import pickle, numpy as np

# Read in chemical evolution tables
global tables, lists
tables = {}
lists = {}
try:
    tables['agb'] = pickle.load(open('agb_enrich.pck'))
except:
    print "Couldn't open agb_enrich.pck"
    pass
lists['agb'] = {}
lists['agb']['iAgeSort'] = np.argsort(tables['agb']['times'][0])
lists['agb']['ageIndices'] = np.arange(len(tables['agb']['times'][0]))
lists['agb']['iZsort'] = np.argsort(tables['agb']['Zs'])
lists['agb']['Zindices'] = np.arange(len(tables['agb']['Zs']))

global snii_dat, isnAgeSort, SNageIndices, isnZsort, SNZindices
try:
    tables['snii'] = pickle.load(open('snii_enrich.pck'))
except:
    print "Couldn't open snii_enrich.pck"
    pass
lists['snii'] = {}
lists['snii']['iAgeSort'] = np.argsort(tables['snii']['times'][0])
lists['snii']['ageIndices'] = np.arange(len(tables['snii']['times'][0]))
lists['snii']['iZsort'] = np.argsort(tables['snii']['Zs'])
lists['snii']['Zindices'] = np.arange(len(tables['snii']['Zs']))

global star_type
star_type = np.dtype({'names':['tform','init_mass','mass','Z'],
                      'formats':['f','f','f','f']})
