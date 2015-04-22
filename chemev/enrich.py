import numpy as np, pickle, pdb
from . import starlifetime as slt, imf
import pdb

mass_weighted_yields = {}
global tables
tables = {}

def stellar_model_to_enrich_time(infile='ww95/lindner99.pck',Zsun_factor=0.02,
                                 outfile='snii_enrich.pck',imf=imf.Chabrier()):

    dat = pickle.load(open(infile))

    el_yield=dat['element_yield']
    ms=dat['masses']
    Zs=dat['Zs']*Zsun_factor
    
    imfTotalMass = imf.cum_mass(0)
    imf_at_ejms = imf.imf(ms)

    yield_rates = {}
    yield_rates['Z']=np.zeros((len(Zs),len(ms)-1))
    times = np.zeros((len(Zs),len(ms)-1))
    for iZ,met in enumerate(Zs):
        mass_times = slt.lifetime(ms,met)
        mass_ranges = ms[1:]-ms[:-1]
        time_ranges = np.vstack((mass_times[1:],mass_times[:-1]))
        times[iZ,:] = np.mean(time_ranges,axis=0)
        time_bin_sizes = np.diff(time_ranges,axis=0)
        Zmass = np.zeros(len(mass_ranges))
        for el in el_yield.keys():
            summed_yields = np.zeros(len(mass_ranges))
            integrands = el_yield[el][iZ,:]*imf_at_ejms
            for ii,inte in enumerate(integrands[:-1]):
                summed_yields[ii] = (inte*mass_ranges[ii]+
                                     integrands[ii+1]*mass_ranges[ii])/2
                if ((el != 'H') and (el !='He') and (el != 'm_ej')):
                    Zmass[ii]+=summed_yields[ii]
            try:
                yield_rates[el][iZ,:] = summed_yields / imfTotalMass / time_bin_sizes
            except KeyError:
                yield_rates[el]=np.zeros((len(Zs),len(ms)-1))
                yield_rates[el][iZ,:] = summed_yields / imfTotalMass / time_bin_sizes
        yield_rates['Z'][iZ,:] = Zmass / imfTotalMass / time_bin_sizes

    pickle.dump({'yield_rates':yield_rates,'times':np.array(times),
                 'time_bin_sizes':time_bin_sizes,'Zs':np.array(Zs)},
                open(outfile,'w'))

    pdb.set_trace()



def make_table(type,time_steps,infile='ww95/lindner99.pck',Zsun_factor=0.02,
               imf=imf.Chabrier()):
    
    sim_time_ranges = np.vstack((time_steps[:-1],time_steps[1:]))
    sim_time_bin_size = np.mean(np.diff(sim_time_ranges,axis=0))
    
    dat = pickle.load(open(infile))

    el_yield=dat['element_yield']
    ms=dat['masses']
    Zs=np.array(dat['Zs'])*Zsun_factor
    isortZs = np.argsort(Zs)
    
    imfTotalMass = imf.cum_mass(0)
    imf_at_ejms = imf.imf(ms)

    yield_rates = {}
    maxoutlen=0
    for iZ,met in enumerate(Zs):
        mass_times = slt.lifetime(ms,met)
        out = np.digitize(mass_times,bins=time_steps)
        maxoutlen = np.max([len(np.unique(out)),maxoutlen])

    maxoutlen = np.min([maxoutlen,len(ms)-1])
    yield_rates['Z']=np.zeros((len(Zs),maxoutlen))
    times = np.zeros((len(Zs),maxoutlen))

    for iZ,met in enumerate(Zs):
        mass_times = slt.lifetime(ms,met)
        out = np.digitize(mass_times,bins=time_steps)
        mass_ranges = ms[1:]-ms[:-1]
        time_ranges = np.vstack((mass_times[1:],mass_times[:-1]))
        time_bin_sizes = np.diff(time_ranges,axis=0)
        Zmass = np.zeros(maxoutlen)
        for el in el_yield.keys():
            summed_yields = np.zeros(maxoutlen)
            integrands = el_yield[el][iZ,:]*imf_at_ejms
            for ii,inte in enumerate(integrands[:-1]):
                if (time_bin_sizes[0,ii] < sim_time_bin_size): itime = out[ii]
                else: itime = maxoutlen-ii
                summed_yields[itime-1] += (inte*mass_ranges[ii]+
                                     integrands[ii+1]*mass_ranges[ii])/2
                if ((el != 'H') and (el !='He') and (el != 'm_ej')):
                    Zmass[itime-1]+=summed_yields[itime-1]

            for it in np.arange(maxoutlen):
                t_range = np.max([time_bin_sizes[0,maxoutlen-it-1], 
                                  sim_time_bin_size])
                try:
                    yield_rates[el][isortZs[iZ],it] = \
                        summed_yields[it] / imfTotalMass / t_range
                except KeyError:
                    yield_rates[el]=np.zeros((len(Zs),maxoutlen))
                    yield_rates[el][isortZs[iZ],it] = \
                        summed_yields[it] / imfTotalMass / t_range
                times[isortZs[iZ],it] = \
                    np.max([time_steps[it+1],mass_times[maxoutlen-it-1]])
                yield_rates['Z'][isortZs[iZ],it]=Zmass[it] /imfTotalMass / t_range

    try:
        tables[type]['yield_rates'] = yield_rates
    except KeyError:
        tables[type] = {}
        tables[type]['yield_rates'] = yield_rates

    tables[type]['times'] = times
    tables[type]['Zs'] = Zs[isortZs]

