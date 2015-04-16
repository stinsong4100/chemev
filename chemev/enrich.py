import numpy as np, pickle, pdb
from . import zones,starlifetime as slt, imf
import pdb

mass_weighted_yields = {}
def stellar_model_to_enrich_time(infile='ww95/lindner99.pck',Zsun_factor=0.02,
                                 outfile='snii_enrich.pck',imf=imf.Chabrier()):

#Limongi only has Zsol and Z=0
#dat = pickle.load(open('limongi12/limongi12.dat'))
    dat = pickle.load(open(infile))

    el_yield=dat['element_yield']
    ms=dat['masses']
    Zs=dat['Zs']*Zsun_factor
    
    imfTotalMass = imf.cum_mass(0)
    imf_at_ejms = imf.imf(ms)

    yield_rates = {}
    times = np.zeros((len(Zs),len(ms)-1))
    for iZ,met in enumerate(Zs):
        mass_times = slt.lifetime(ms,met)
        mass_ranges = ms[1:]-ms[:-1]
        time_ranges = np.vstack((mass_times[1:],mass_times[:-1]))
        times[iZ,:] = np.mean(time_ranges,axis=0)
        time_bin_sizes = np.diff(time_ranges,axis=0)
        for el in el_yield.keys():
            summed_yields = np.zeros(len(mass_ranges))
            integrands = el_yield[el][iZ,:]*imf_at_ejms
            for ii,inte in enumerate(integrands[:-1]):
                summed_yields[ii] = (inte*mass_ranges[ii]+
                                     integrands[ii+1]*mass_ranges[ii])/2
            try:
                yield_rates[el][iZ,:] = summed_yields / imfTotalMass / time_bin_sizes
            except KeyError:
                yield_rates[el]=np.zeros((len(Zs),len(ms)-1))
                yield_rates[el][iZ,:] = summed_yields / imfTotalMass / time_bin_sizes

    pickle.dump({'yield_rates':yield_rates,'times':np.array(times),'Zs':np.array(Zs)},
                open(outfile,'w'))

    pdb.set_trace()

