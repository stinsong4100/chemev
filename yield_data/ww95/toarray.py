import numpy as np, pickle, pdb

infiles = ['Z0.txt','Z1e-4sol.txt','Z0.01.txt','Z0.1sol.txt','Zsol.txt']
Zs = [0,2e-6,2e-4,2e-3,0.02]
mass_names = ('11A','12A','13A','15A','18A','19A','20A','22A','25A','25B','30A','30B','35A','35B','35C','40A','40B','40C')
formats= (len(mass_names)*['f'])
mass_type = np.dtype({'names':mass_names,'formats':formats})

# load decay products
pre,post = np.genfromtxt('../general/decay_products.txt',skip_header=1,dtype='string',unpack=True)
decay_products = {}
for ip,pre in enumerate(pre): decay_products[pre]=post[ip]

isotope_yield={}
element_yield={}

el_trans = {'1H':'^1^H','2H':'^2^H','3He':'^3^He','4He':'^4^He','7Li':'^7^Li',
            '7Be':'^7^Be','9Be':'^9^Be','m_ej':'m_ej'}
isos = np.genfromtxt('isonames.txt',usecols=(0),dtype='string')

for iz,infile in enumerate(infiles):
    
    prods = np.genfromtxt(infile,names=True)
    
    for ii,iso in enumerate(isos):
        if iso in el_trans.keys(): ciso = el_trans[iso]
        else:  ciso = '^'+iso[:2]+'^'+iso[2:]

        if ciso not in decay_products: 
            continue
        if decay_products[ciso] not in isotope_yield: 
            isotope_yield[decay_products[ciso]] = np.zeros((len(Zs)),dtype=mass_type)
        for m in prods.dtype.fields.keys():
            isotope_yield[decay_products[ciso]][m][iz] += prods[m][ii]
    if 'm_ej' not in isotope_yield: 
        isotope_yield['m_ej'] = np.zeros((len(Zs)),dtype=mass_type)
    for m in prods.dtype.fields.keys():
        isotope_yield['m_ej'][m][iz] += prods[m][ii-1]
        
    for iso in isotope_yield.keys():
        element = iso.split('^')[-1]
        if iz==0: element_yield[element] = np.zeros((len(Zs)),dtype=mass_type)
        for m in mass_names:
            element_yield[element][m][iz] += isotope_yield[iso][m][iz]
    if 'm_ej' not in element_yield: 
        element_yield['m_ej'] = np.zeros((len(Zs)),dtype=mass_type)
    for m in prods.dtype.fields.keys():
        element_yield['m_ej'][m][iz] += prods[m][ii-1]

    
pdb.set_trace()
pickle.dump({'isotope_yield':isotope_yield,'element_yield':element_yield,
             'masses':mass_names,'Zs':np.array(Zs)},open('ww95.pck','w'))
