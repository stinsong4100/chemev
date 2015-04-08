import numpy as np, pickle, pdb
from chemev import starlifetime as slt
import chemev

mass_weighted_yields = {}
imf = chemev.imf.Chabrier()

dat = pickle.load(open('limongi12/limongi12.dat'))

el_yield=dat['element_yield']
masses=dat['masses']
Zs=dat['Zs']
iZ=1

output_times = np.arange(slt.lifetime(80,Zs[iZ]*0.02),slt.lifetime(8,Zs[iZ]*0.02),1e6)
time_masses=slt.starMass(output_times,0.02)
imfTotalMass = imf.cum_mass(0)
cum_masses = imf.cum_mass(time_masses)
mass_fractions = (cum_masses[1:] - cum_masses[:-1]) / imfTotalMass
imf_at_masses = imf.imf(time_masses)

for el in el_yield.keys():
    yield_at_masstimes = np.interp(time_masses,masses,el_yield[el][iZ,:]/masses)
    summed_yields = np.trapz(yield_at_masstimes*imf_at_masses,x=np.vstack([time_masses[1:],
                                                            time_masses[:-1]]).T)
    mass_weighted_yields[el] = summed_yields
    pdb.set_trace()

pickle.dump({'mass_weighted_yields':mass_weighted_yields,'times':output_times},
            open('time_yields.pck','w'))
