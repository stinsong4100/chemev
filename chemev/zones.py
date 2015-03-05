"""
zones
====

Implements the zone class, which is the fundamental container for gas
in the chemical evolution models.
"""

import numpy as np
import logging
from . import star.Star

class Zone():

    """
    Generic class representing a zone.
    """

    def __init__(self,radius,mass,Z,abunds,*args):
        self.r = radius
        self.mass = mass
        self.abunds = abunds
        self.Z = Z
        self._stars = np.array([])

    def add_star():
        """
        Adds a new star to the zone.
        """
        
        self._stars.append(Star(mass,now,self.Z,self.abunds))

    def enrich():
        """
        Loop through all the stars in the zone and enrich all the elements
        we are tracing for one timestep.
        """
        self.abunds,self.Z += snii(tnow, self._stars)
        self.abunds,self.Z += snia(tnow, self._stars)
        self.abunds,self.Z += agb(tnow, self._stars)

    def outflow():
        """
        Move material out of this zone to target zone
        """

    def inflow():
        """
        Add mass from another zone.
        """


def create_zones(n):
    """
    Initializes zone creation, automates creating zones in a disk.
    """
    
    disk_zones = [Zone(i*dr) for i in range(n)]
