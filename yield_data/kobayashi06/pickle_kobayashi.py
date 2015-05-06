import numpy as np, pickle, pdb

infiles = ['Z0.txt','Z0.001.txt','Z0.004.txt','Z0.02.txt']
Zs = [0 ,0.001,0.004,0.02]
masses = np.array([13,15,18,20,25,30,40])

# load decay products
pre,post = np.genfromtxt('../general/decay_products.txt',skip_header=1,dtype='string',unpack=True)
decay_products = {}
for ip,pre in enumerate(pre): decay_products[pre]=post[ip]

isotope_yield={}
element_yield={}

isos = np.genfromtxt('isonames.txt',usecols=(0),dtype='string')

for iz,infile in enumerate(infiles):
    
    prods = np.genfromtxt(infile)
    
    for ii,iso in enumerate(isos):
        if iso not in decay_products: continue
        if decay_products[iso] not in isotope_yield: 
            isotope_yield[decay_products[iso]] = np.zeros((len(Zs),len(masses)))
        isotope_yield[decay_products[iso]][iz,:] += prods[ii]

    if 'm_ej' not in isotope_yield: 
        isotope_yield['m_ej'] = np.zeros((len(Zs),len(masses)))
    isotope_yield['m_ej'][iz,:] = masses - prods[1]
        
for iso in isotope_yield.keys():
    if iso =='m_ej': continue
    element = iso.split('^')[-1]
    if element not in element_yield: element_yield[element] = np.zeros((len(Zs),len(masses)))
    element_yield[element] += isotope_yield[iso]
    
element_yield['m_ej'] = isotope_yield['m_ej']

    
pickle.dump({'isotope_yield':isotope_yield,'element_yield':element_yield,
             'masses':masses,'Zs':np.array(Zs)},open('ww95.pck','w'))

pdb.set_trace()
