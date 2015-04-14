import numpy as np, pickle, pdb

infiles = ['Z_0.0001.txt','Z_0.004.txt','Z_0.008.txt','Z_0.02.txt']
Zs = [1e-4,4e-3,8e-3,0.02]
masses = np.array([ 1.  ,  1.25,  1.5 ,  1.75,  1.9 ,  2.  ,  2.25,  2.5 ,  3.  ,
                    3.5 ,  4.  ,  4.5 ,  5.  ,  5.5 ,  6.  ,  6.5 ])

# load decay products
pre,post = np.genfromtxt('../general/decay_products.txt',skip_header=1,dtype='string',unpack=True)
decay_products = {}
for ip,pre in enumerate(pre): decay_products[pre]=post[ip]

isotope_yield={}
element_yield={}

el_trans = {'p':'^1^H','d':'^2^H','he3':'^3^He','he4':'^4^He','li7':'^7^Li',
            'be7':'^7^Be','b8':'^8^B'}

for iz,infile in enumerate(infiles):
    
    mis,mets,mfs,prods = np.genfromtxt(infile,usecols=(0,1,2,6),unpack=True)
    isos = np.genfromtxt(infile,usecols=(3),dtype='string')
    dum_masses, imasses = np.unique(mis,return_index=True)
    Z = np.unique(mets)
    iZ = np.argwhere(Zs==Z)[0]

    for ii,iso in enumerate(isos):
        if iso in el_trans.keys(): ciso = el_trans[iso]
        else: 
            capiso = iso.capitalize()
            ciso = '^'+capiso[-2:]+'^'+capiso[:-2]

        if ciso not in decay_products: continue
        if decay_products[ciso] not in isotope_yield: 
            isotope_yield[decay_products[ciso]] = np.zeros((len(Zs),len(masses)))
        try:  iMass = np.argwhere(masses == mis[ii])[0]
        except IndexError:
            if mis[ii] == 2.1: iMass = np.argwhere(masses==2)[0]
        isotope_yield[decay_products[ciso]][iZ,iMass] += prods[ii]
    if 'm_ej' not in isotope_yield: 
        isotope_yield['m_ej'] = np.zeros((len(Zs),len(masses)))
    isotope_yield['m_ej'][iZ,:len(dum_masses)] = dum_masses - mfs[imasses]
        
    for iso in isotope_yield.keys():
        element = iso.split('^')[-1]
        if iz==0: element_yield[element] = np.zeros((len(Zs),len(masses)))
        element_yield[element] += isotope_yield[iso][iZ,:]

element_yield['m_ej'] = isotope_yield['m_ej']

pickle.dump({'isotope_yield':isotope_yield,'element_yield':element_yield,
             'masses':masses,'Zs':Zs},open('karakas2010.pck','w'))
