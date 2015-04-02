import numpy as np, pickle, pdb

infiles = ['nonrotateZ0.txt','nonrotateZsol.txt']
Zs = [0,1]
masses = [13,15,20,25,30,35,50,80]

# load decay products
pre,post = np.genfromtxt('decay.txt',skip_header=1,dtype='string',unpack=True)
decay_products = {}
for ip,pre in enumerate(pre): decay_products[pre]=post[ip]

isotope_yield={}
element_yield={}

for iz,infile in enumerate(infiles):
    
    prods = np.genfromtxt(infile,usecols=np.arange(1,len(masses)+1))
#,names=['%d m'%m for m in masses])
    isos = np.genfromtxt(infile,usecols=(0),dtype='string')
    
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
