import numpy as np, pickle, pdb

infiles = ['Z0.txt','Z1e-4sol.txt','Z0.01.txt','Z0.1sol.txt','Zsol.txt']
Zs = [0,2e-6,2e-4,2e-3,0.02]
#mass_names = ('11A','12A','13A','15A','18A','19A','20A','22A','25A','25B','30A','30B','35A','35B','35C','40A','40B','40C')
#formats= (len(mass_names)*['f'])
#mass_type = np.dtype({'names':mass_names,'formats':formats})
masses = np.array([11,12,13,15,18,19,20,22,25,30,35,40])
models = ['A','A','A','A','A','A','A','A','A','B','C','C']

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
            isotope_yield[decay_products[ciso]] = np.zeros((len(Zs),len(masses)))
        mass_bools = np.zeros(len(masses),dtype='bool')
        for m in prods.dtype.fields.keys():
            imass = np.argwhere(int(m[:-1])==masses)[0][0]
            # take highest energy model available
            if m[-1] == models[imass]:
                isotope_yield[decay_products[ciso]][iz,imass] += prods[m][ii]
                mass_bools[imass]=True
        # interpolate / extrapolate fields not filled
        igood = np.argwhere(mass_bools)
        ibad = np.argwhere(mass_bools==False)
        if len(ibad) == 0: continue
        for intm in ibad.T[0]:
            if masses[intm] > 12:  #interpolate
                isotope_yield[decay_products[ciso]][iz,intm] = \
                    np.interp(masses[intm], masses[igood].T[0], 
                              isotope_yield[decay_products[ciso]][iz,igood].T[0])
            else: #extrapolate
                fit = np.polyfit(masses[igood].T[0], 
                                 isotope_yield[decay_products[ciso]][iz,igood].T[0],2)
                line = np.poly1d(fit)
                isotope_yield[decay_products[ciso]][iz,intm] = \
                    np.max([0,line(masses[intm])])

    if 'm_ej' not in isotope_yield: 
        isotope_yield['m_ej'] = np.zeros((len(Zs),len(masses)))
    for m in prods.dtype.fields.keys():
        if m[-1] == 'A':
            imass = np.argwhere(int(m[:-1])==masses)[0][0]
            isotope_yield['m_ej'][iz,imass] = prods[m][ii]
        
for iso in isotope_yield.keys():
    if iso =='m_ej': continue
    element = iso.split('^')[-1]
    if element not in element_yield: element_yield[element] = np.zeros((len(Zs),len(masses)))
    element_yield[element] += isotope_yield[iso]
    
element_yield['m_ej'] = isotope_yield['m_ej']

    
pickle.dump({'isotope_yield':isotope_yield,'element_yield':element_yield,
             'masses':masses,'Zs':np.array(Zs)},open('ww95.pck','w'))

pdb.set_trace()
