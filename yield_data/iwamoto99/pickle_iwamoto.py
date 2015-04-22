import numpy as np, pickle, pdb

isotope_yield={}
element_yield={}
infile='iwamoto99.txt'

prods = np.genfromtxt(infile,usecols=(2),skip_header=2)
isos = np.genfromtxt(infile,usecols=(0),dtype='string',skip_header=2)

for ii,iso in enumerate(isos):
    ciso = '^'+iso[:2]+'^'+iso[2:]
    isotope_yield[ciso] = prods[ii]
        
for iso in isotope_yield.keys():
    element = iso.split('^')[-1]
    try:
        element_yield[element] += isotope_yield[iso]
    except:
        element_yield[element] = isotope_yield[iso]

pickle.dump({'element_yield':element_yield,'isotope_yield':isotope_yield},
            open('iwamoto99sniaW7.pck','w'))
