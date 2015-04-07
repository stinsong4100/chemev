import numpy as np, pickle, pdb
from chemev import starlifetime as slt

dat = pickle.load(open('limongi12/limongi12.dat'))

el_yield=dat['element_yield']
masses=dat['masses']
Zs=dat['Zs']
iZ=1

output_times = np.arange(slt.lifetime(80,iZ*0.02),slt.lifetime(8,iZ*0.02),1e6)
time_masses=slt.starMass(output_times,0.02)

for el in el_yield.keys():
    yield_at_masstimes = np.interp(time_masses,masses,el_yield[el][iZ,:])
    summed_yields = np.trapz(yield_at_masstimes,x=np.vstack([time_masses[1:],
                                                            time_masses[:-1]]).T)

pdb.set_trace()
