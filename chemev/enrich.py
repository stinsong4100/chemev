import numpy as np, pickle, pdb
from . import starlifetime as slt, imf, snia
import pdb

mass_weighted_yields = {}
global tables
tables = {}
min_SNII_mass = 8.0
dMetals = 0.02

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

#    pdb.set_trace()



def make_table(type,time_steps,infile='yields/snii/kobayashi.pck',
               imf=imf.Chabrier(),plot=False):
    
    sim_time_ranges = np.vstack((time_steps[:-1],time_steps[1:]))
    sim_time_bin_size = np.mean(np.diff(sim_time_ranges,axis=0))

    added_min_mass=False
    
    dat = pickle.load(open(infile))

    el_yield=dat['element_yield']
    ms=dat['masses']
        
    Zs=np.array(dat['Zs'])   #Zs in mass fraction, NOT fraction of solar
    isortZs = np.argsort(Zs)
    
    if plot:
        import pylab as plt
        colors=['r','g','b','m','k','y','c']
        nplots = len(el_yield)
        ncols = 6
        nrows = np.int(np.ceil(np.float(nplots) / ncols))
        f,ax = plt.subplots(nrows,ncols,sharex=True,figsize=(12,nrows*2+2))
        f.subplots_adjust(top=0.99,hspace=0,bottom=0.04,right=0.99)

        m_range = np.arange(min_SNII_mass,50,1)

    if (((type == 'snii') and (np.min(ms) > min_SNII_mass)) or plot):
        extr_ms = np.arange(min_SNII_mass,np.min(ms),1)
        for iel, el in enumerate(el_yield.keys()):
            extrap_yield = np.zeros((len(Zs),len(extr_ms)))
            new_yields = np.zeros((len(Zs),len(ms)+len(extr_ms)))

            for iZ,met in enumerate(Zs):
                fit = np.polyfit(ms,el_yield[el][iZ,:],2)
                line = np.poly1d(fit)
                new_yields[iZ,:] = np.concatenate([np.maximum(0,line(extr_ms)),el_yield[el][iZ,:]])
                if plot: 
                    c = colors[iZ]
                    ir = np.int(np.floor(iel/ncols))
                    ic = np.mod(iel,ncols)
                    ax[ir,ic].plot(ms,el_yield[el][iZ,:],'o',color=c)
                    ax[ir,ic].text(0.02,0.85,el,transform=ax[ir,ic].transAxes)
                    ax[ir,ic].plot(m_range,line(m_range),color=c)
                    ax[ir,ic].plot(extr_ms,line(extr_ms),'o',color=None,mfc=None,mec=c)
            el_yield[el] = new_yields
        ms = np.insert(ms,0,extr_ms)
#        pdb.set_trace()

    if plot: f.savefig('elements.png')

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

#    pdb.set_trace()
    tables[type]['times'] = times
    tables[type]['Zs'] = Zs[isortZs]

def make_snia_table(time_steps,snia_model='Maoz',
                    infile='yields/snia/iwamoto.pck',
                    imf=imf.Chabrier(),min_time_step=1e8):

    dat = pickle.load(open(infile))
    el_yield=dat['element_yield']
    Zs  = np.array(dat['Zs'])
    isortZs = np.argsort(Zs)

    maxSNII_lifetime = slt.lifetime(min_SNII_mass,dMetals)

    if min_time_step == 'log':
        snia_ts = np.logspace(np.log10(maxSNII_lifetime), 
                              np.log10(np.max(time_steps)),50)
        time_ranges = np.vstack((snia_ts[:-1],snia_ts[1:]))
        time_step = np.diff(time_ranges,axis=0)
    else:
        sim_time_ranges = np.vstack((time_steps[:-1],time_steps[1:]))
        sim_time_bin_size = np.mean(np.diff(sim_time_ranges,axis=0))
        time_step = np.max([sim_time_bin_size,min_time_step])
        snia_ts = np.arange(maxSNII_lifetime, np.max(time_steps),time_step)
    
    dMtot = imf.cum_mass(0)

    nsnias = np.zeros(len(snia_ts))

    print 'Calculating %d SNIa rates'%(len(snia_ts))
    if snia_model == 'Greggio1993':
        for it,time in enumerate(snia_ts):
            dMStarMin = slt.starMass(time+time_step, dMetals)
            dMStarMax = slt.starMass(time, dMetals)
            
            nsnias[it] = snia.NSNIa(dMStarMin,dMStarMax,imf) / dMtot

    elif snia_model == 'Maoz':
        # Definite integral of 3.75e-3*ts**-1.12 for all time intervals
        nsnias=-3.75e-3/0.12*snia_ts[1:]**-0.12 + \
            3.75e-3/0.12*snia_ts[:-1]**-0.12

    yield_rates = {}
    for el in dat['element_yield'].keys():
        yield_rates[el] = nsnias*dat['element_yield'][el][np.newaxis].T / time_step

    tables['snia'] = {}
    tables['snia']['yield_rates'] = yield_rates
    tables['snia']['times'] = snia_ts
    tables['snia']['Zs'] = Zs[isortZs]
