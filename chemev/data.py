import pickle, numpy as np

# Read in chemical evolution tables
try:
    global agb_dat 
    agb_dat = pickle.load(open('agb_enrich.pck'))
    global snii_dat 
    snii_dat = pickle.load(open('snii_enrich.pck'))
except:
    pass

global star_type
star_type = np.dtype({'names':['tform','init_mass','mass','Z'],
                      'formats':['f','f','f','f']})
