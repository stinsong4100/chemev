import numpy as np, pickle, pdb

infiles = ['Z_0.0001.txt','Z_0.004.txt','Z_0.008.txt','Z_0.02.txt']
Zs = [1e-4,4e-3,8e-3,0.02]

# load decay products
pre,post = np.genfromtxt('../general/decay_products.txt',skip_header=1,dtype='string',unpack=True)
decay_products = {}
for ip,pre in enumerate(pre): decay_products[pre]=post[ip]

isotope_yield={}
element_yield={}

el_trans = {'P':'^1^H','D':'^2^H','He3':'^3^He','He4':'^4^He','Li7':'^7^Li',
            'Be7':'^7^Be','B8':'^8^B'}

for iz,infile in enumerate(infiles):
    
    prods = np.genfromtxt(infile,usecols=np.arange(1,len(masses)+1))
    isos = np.genfromtxt(infile,usecols=(3),dtype='string').capitalize()
    masses = np.unique(np.genfromtxt(infile,usecols=(0)))
    
    for ii,iso in enumerate(isos):
        if iso not in decay_products: continue
        if decay_products[iso] not in isotope_yield: 
            isotope_yield[decay_products[iso]] = np.zeros((len(Zs),len(masses)))
        isotope_yield[decay_products[iso]][iz,:] += prods[ii]
        
    for iso in isotope_yield.keys():
        element = iso.split('^')[-1]
        if iz==0: element_yield[element] = np.zeros((len(Zs),len(masses)))
        element_yield[element] += isotope_yield[iso][iz,:]

pickle.dump({'isotope_yield':isotope_yield,'element_yield':element_yield,
             'masses':masses,'Zs':Zs},open('limongi12.dat','w'))
