"""
star
====

Implements the star class, the star data container object
"""

import numpy as np
import logging

class Star():

    """
    Defines star class.
    """

    def __init__(self,tform,Z,abunds,sf_mode='sim',mass=None,zonemass=None,zonearea=None):
        self.tform = tform
        self.Z = Z
        self.abunds = abunds

        if sf_mode=='sim': self.mass=mass
        elif sf_mode=='kennicutt': self.mass = star_form_efficiency * zonemass #/zonearea

        self.init_mass = self.mass

            
