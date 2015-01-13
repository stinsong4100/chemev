#!/usr/bin/env python

"""
FILE
    plot_ofe_mdf.py

DESCRIPTION
    Plot [O/Fe]--[Fe/H], MDF, and [O/Fe] distribution function.

"""

version = 0.0

import os
import sys
home = os.path.expanduser('~')
stem_flexce = home + '/FlexCE_v' + str(version) + '/'
stem_scripts = stem_flexce + 'scripts/'
sys.path.append(stem_scripts)

import numpy as np
import pickle
import traceback
import collections
import matplotlib.pyplot as plt

import flexce_main

import plotting_scripts
from plotting_scripts import PlotSim


##### Set path #####
stem_runs = stem_flexce + 'runs/'
stem_param = stem_flexce + 'param_files/'
stem_plot = stem_flexce + 'plots/'
####################


##### Load Chemical Evolution Models #####
for i in [0]:
    try:
        for ftype in ['box', 'ab']:
            file_in = stem_runs + ftype + '%s.pck' % i
            fin = open(file_in, 'r')
            globals()[ftype + str(i)] = pickle.load(fin)
            fin.close()
    except (KeyError, EOFError) as e:
        print 'sim' + str(i) + ' did not load'

##########################################

names = [
    r'box0: fiducial',
    ]
sim_ids = [name.split(':')[0].lstrip('box') for name in names]
labels = [name.split(': ')[1] for name in names]

colors = ['b', 'orange', 'm', 'r']

s1_tmp = {}
for sim_id, label, color in zip(sim_ids, labels, colors):
    try:
        s1_tmp['sim' + sim_id] = {'box': globals()['box' + sim_id],            
        'ab': globals()['ab' + sim_id], 'name': label, 'c': color}
    except KeyError as e:
        print traceback.print_exc()

suite1 = collections.OrderedDict(sorted(s1_tmp.items()))

ps1 = PlotSim(suite1)
ps1.plot_xfe_mdf('O', runs=['sim0'],
                 xylim=[-2, 0.49, -0.1, 0.5],
                 mdf_xylim=[-2, 0.5, 0, 10],
                 xfe_xylim=[0, 20, -0.1, 0.5],
                 xfe_args={'min': -0.1, 'max': 0.5, 'delta': 0.001}, 
                 time_pts=True, tlab=True, time_pts_offset=[0.03, 0.02])
plt.savefig(stem_plot + 'ofe_mdf.pdf')
