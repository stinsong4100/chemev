"""
snii
====

Contains all the methods that do SNII enrichment.
"""

import numpy as np
import logging

def snii(tnow, stars):
    print 'do nothing'

def ira(tnow, stars):
    """
    Instantaneous recycling approximation
    """
    istars = np.where(((tnow - stars.tform)>3.5e6) and 
                      ((tnow - stars.tform)<3.5e6+tdelta))

    
