import numpy as np, pickle, pdb

infiles = ['Z0.txt','Z1e-4sol.txt','Z0.01.txt','Z0.1.txt','Zsol.txt']
Zs = [0,2e-6,2e-4,2e-3,0.02]
masses = [11,12,13,15,18,19,20,22,25,30,30,35,35,35,40,40]

# load decay products
pre,post = np.genfromtxt('../general/decay_products.txt',skip_header=1,dtype='string',unpack=True)
decay_products = {}
for ip,pre in enumerate(pre): decay_products[pre]=post[ip]

isotope_yield={}
element_yield={}

el_trans = {'1H':'^1^H','2H':'^2^H','3He':'^3^He','4He':'^4^He','7Li':'^7^Li',
            '7Be':'^7^Be','9Be':'^9^Be'}

for iz,infile in enumerate(infiles):
    
    prods = np.genfromtxt(infile,usecols=np.arange(1,len(masses)+1),names=True)
#,names=['%d m'%m for m in masses])
    isos = np.genfromtxt(infile,usecols=(0),dtype='string')
    
    for ii,iso in enumerate(isos):
        if iso in el_trans.keys(): ciso = el_trans[iso]
        else: 
            ciso = '^'+iso[-2:]+'^'+iso[:-2]

        if ciso not in decay_products: continue
        if decay_products[ciso] not in isotope_yield: 
            isotope_yield[decay_products[ciso]] = np.zeros((len(Zs),len(masses)))
        isotope_yield[decay_products[ciso]][iz,:] += prods[ii]
        
    for iso in isotope_yield.keys():
        element = iso.split('^')[-1]
        if iz==0: element_yield[element] = np.zeros((len(Zs),len(masses)))
        element_yield[element] += isotope_yield[iso][iz,:]

pickle.dump({'isotope_yield':isotope_yield,'element_yield':element_yield,
             'masses':masses,'Zs':Zs},open('ww95.pck','w'))
