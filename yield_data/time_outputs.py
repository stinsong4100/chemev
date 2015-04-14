import numpy as np, pickle, pdb
from chemev import starlifetime as slt
import chemev

mass_weighted_yields = {}
imf = chemev.imf.Chabrier()

#Limongi only has Zsol and Z=0
#snii_dat = pickle.load(open('limongi12/limongi12.dat'))
snii_dat = pickle.load(open('ww95/lindner99.dat'))
agb_dat = pickle.load(open('karakas10/karakas2010.pck'))

snii_el_yield=snii_dat['element_yield']
snii_ms=snii_dat['masses']
snii_Zs=snii_dat['Zs']*0.02

agb_el_yield=agb_dat['element_yield']
agb_ms=agb_dat['masses']
agb_Zs=agb_dat['Zs']

ej_ms=np.concatenate((agb_ms,snii_ms))

#output_times = np.arange(slt.lifetime(80,Zs[iZ]*0.02),slt.lifetime(8,Zs[iZ]*0.02),1e6)
#time_masses=slt.starMass(output_times,0.02)

imfTotalMass = imf.cum_mass(0)
imf_at_ejms = imf.imf(ej_ms)

yield_rates = {}
good_Zs = []
nZs = np.min([len(snii_Zs),len(agb_Zs)])
times = np.zeros((nZs,len(ej_ms)-1))
for iagbZ,met in enumerate(agb_Zs):
    try:
        isniiZ = np.argwhere(snii_Zs == met)[0][0]
    except:
        continue
    good_Zs.append(met)
    mass_times = slt.lifetime(ej_ms,met)
    time_ranges = np.vstack((mass_times[1:],mass_times[:-1]))
    times[iagbZ,:] = np.mean(time_ranges,axis=0)
    time_bin_sizes = np.diff(time_ranges,axis=0)
    agb_els = np.array(agb_el_yield.keys())
    good_els = agb_els[np.in1d(agb_el_yield.keys(),snii_el_yield.keys())]
    for el in good_els:
        el_yield = np.hstack((agb_el_yield[el][iagbZ,:],snii_el_yield[el][isniiZ,:]))
        summed_yields = np.trapz(el_yield*imf_at_ejms,
                                 x=np.vstack([ej_ms[:-1],ej_ms[1:]]).T)
        try:
            yield_rates[el][iagbZ,:] = summed_yields / imfTotalMass / time_bin_sizes
        except KeyError:
            yield_rates[el]=np.zeros((nZs,len(ej_ms)-1))
            yield_rates[el][iagbZ,:] = summed_yields / imfTotalMass / time_bin_sizes


pickle.dump({'yield_rates':yield_rates,'times':times,'Zs':good_Zs},
            open('time_yields.pck','w'))
