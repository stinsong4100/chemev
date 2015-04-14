import numpy as np, pickle, pdb

elements = ['H','He','C','N','O','Mg','Al','Si','S','Ar','Ca','Fe','Ni','Zn','Cr','Mn']

element_yield={}

prods = np.genfromtxt('lindnerSummary.txt',usecols=np.arange(2,len(elements)+2),
                      names=elements,comments='A',skip_footer=26)
Z,mass = np.genfromtxt('lindnerSummary.txt',usecols=(0,1),unpack=True,comments='A',
                       skip_footer=26)
Zs = np.unique(Z)
nZs = len(Zs)
masses = np.unique(mass)
metals = []
nms = len(masses)
element_yield['m_ej'] = np.zeros((nZs,nms))
for ie,el in enumerate(elements):
    element_yield[el] = np.zeros((nZs,nms))
    for iz in np.arange(nZs):
        if ie==0: metals.append( np.unique(Z[iz*nms:(iz+1)*nms])[0])
        element_yield[el][iz,:] += prods[el][iz*nms:(iz+1)*nms]
    element_yield['m_ej'] += element_yield[el]

pickle.dump({'element_yield':element_yield,
             'masses':masses,'Zs':np.array(metals)},open('lindner99.dat','w'))
