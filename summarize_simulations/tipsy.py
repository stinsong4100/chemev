"""
tipsy
=====

Creates summary of evolution of pkdgrav / ChaNGa / gasoline simulations
"""

import sys, numpy as np, pdb, pynbody, pickle, glob

sim_time=13.73e9
step_time=1e6

# Star formation history from starlog
sl = pynbody.tipsy.StarLog('g1536.starlog')
sfh,startimes = np.histogram(sl['tform'].in_units('Gyr'),weights=sl['massform'].in_units('Msol'),
                             bins=sim_time/step_time)

gastimes=[]
gasmass=[]
# Gas Disk Mass evolution from data files
dfs = glob.glob('/disk2/stinson/g1536/changa/esn1.5/00????/g1536.00????.data')
for df in dfs:
    dat = pickle.load(open(df))
    gastimes.append(dat['time'])
    gasmass.append(dat['mdiskcoolgas'])


pickle.dump({'sfh':sfh,'startimes':startimes,'gastimes':gastimes,'gasmass':gasmass},
            open('g1536_sfh+gas.pck','w'))
