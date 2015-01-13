"""
FILE
    setup.py

DESCRIPTION

    Generate yield grids.
    
"""

version = 0.0

import os
import sys
import numpy as np

home = os.path.expanduser('~')
stem_flexce = home + '/FlexCE_v' + str(version) + '/'
stem_yields = stem_flexce + 'yields/'
stem_calc_yields = stem_flexce + 'scripts/calc_yields/'


make_ww95 = False

if len(sys.argv) > 1:
    if sys.argv[1] == 'ww95':
        make_ww95 = True
    else:
        print ' '.join(['\nRun without keyword or with one of the following',
                       'accepted keywords:\nww95\n'])
        sys.exit(1)



print '\nGenerating pickled yield grids...\n'


ylds = {
    'busso01': {},
    'cescutti06': {},
    'iwamoto99': {},
    'karakas10': {},
    'limongi06': {},
    'sneden08': {},
    'ww95': {}
    }

for k in np.sort(ylds.keys()):
    py_file = k + '_yields.py'
    pck_file = 'interp_yields.pck'
    stem = k + '/'
    if k == 'busso01':
        pck_file = 'busso01_yields.pck'
    elif k == 'cescutti06':
        pck_file = 'cescutti06_yields.pck'
    elif k == 'iwamoto99':
        pck_file = 'w70_yields.pck'
    elif k == 'karakas10':
        stem = stem + 'iso_yields/'
    if k == 'limongi06':
        py_file = 'limongi_chieffi_yields.py'
        stem = stem + 'iso_yields/'
    elif k == 'sneden08':
        py_file = 'sneden08.py'
        pck_file = 'sneden08.pck'
        stem = 'general/'
    elif k == 'ww95':
        stem = stem + 'half_fe/'
    ylds[k]['script'] = py_file
    ylds[k]['pck'] = pck_file
    ylds[k]['stem'] = stem


for k in np.sort(ylds.keys()):
    if os.path.isfile(stem_yields + ylds[k]['stem'] + ylds[k]['pck']):
        print 'yields already exist: ', k
    else:
        try:
            if k == 'ww95':
                if not make_ww95:
                    continue
            print 'generating yields: ', k
            os.system('python ' + stem_calc_yields + ylds[k]['script'])
        except:
            print 'couldn\'t generate yields: ', k


print '\n\n'












