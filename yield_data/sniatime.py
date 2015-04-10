import totradpres, startime, snia, numpy as np, pdb
import pylab as pl


dtime = 867000
nsamp = 13700
dMetals = 0.02
pmass = 5.5e4
times = np.zeros(nsamp)
snias = np.zeros(nsamp)

try:
    times,snias = np.loadtxt('snia.dat',unpack=True)
except:
    f = open('snia.dat','w')
    print "time[yr] N_SNIa"
    for i in np.arange(nsamp):
        time = i*dtime
        times[i] = time
        
    # total mass in stars integrated over IMF 
        dMtot = totradpres.calc_mass(0.1)
        
        if (time > 3e7):
            dMStarMin = startime.dSTMStarLtime(time+dtime, dMetals)
            dMStarMax = startime.dSTMStarLtime(time, dMetals)
            
            snias[i] = snia.NSNIa(dMStarMin,dMStarMax) / dMtot
            if np.mod(i,1000)==0: print time, snias[i]
            f.write('%.4e  %.4e\n'%(time,snias[i]) )
    f.close()

pl.clf()
pl.loglog(times,snias/dtime,label='Greggio & Renzini (1983)')

dtds=np.array([2.2e8,1.3e9,7e9])
sniarate=np.array([1.5e-12,3e-13,2e-14])
pl.plot(dtds,sniarate,'or',ms=15,mec='r',label='Maoz+ (2012)')

ts = np.logspace(7.6,10.1)
mdtd=3.75e-3*ts**-1.12
ts=np.insert(ts,0,10**7.5)
mdtd=np.insert(mdtd,0,1e-16)
pl.plot(ts,mdtd,label=r'3.75$\times10^{-3}$ t$^{-1.12}$')

ots = np.logspace(8,10.1)
odtd=3.75e-3*ots**-1.12
ots=np.insert(ots,0,10**7.99)
odtd=np.insert(odtd,0,1e-16)
pl.plot(ots,odtd,'m',label=r'start at 100 Myr')

#
pdb.set_trace()

pl.xlabel('Time [yr]')
pl.ylabel('N$_{SNIa}$ [M$_{\odot}$$^{-1}$ yr$^{-1}$]')
pl.xlim(1e7,2e10)
pl.legend(loc=0)
pl.savefig('sniadtd.png')
