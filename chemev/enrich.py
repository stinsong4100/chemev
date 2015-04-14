import numpy as np, pickle, pdb
from . import starlifetime as slt, imf

mass_weighted_yields = {}
def stellar_model_to_enrich_time(infile='ww95/linder99.dat',Zsun_factor=0.02
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
        time_ranges = np.vstack((mass_times[1:],mass_times[:-1]))
        times[iZ,:] = np.mean(time_ranges,axis=0)
        time_bin_sizes = np.diff(time_ranges,axis=0)
        for el in el_yield.keys():
            summed_yields = np.trapz(el_yield[el][iZ,:]*imf_at_ejms,
                                     x=np.vstack([ms[:-1],ms[1:]]).T)
            try:
                yield_rates[el][iZ,:] = summed_yields / imfTotalMass / time_bin_sizes
            except KeyError:
                yield_rates[el]=np.zeros((len(Zs),len(ms)-1))
                yield_rates[el][iZ,:] = summed_yields / imfTotalMass / time_bin_sizes

    pickle.dump({'yield_rates':yield_rates,'times':times,'Zs':Zs},
            open(outfile,'w'))
