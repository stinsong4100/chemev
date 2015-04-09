"""
zones
====

Implements the zone class, which is the fundamental container for gas
in the chemical evolution models.
"""

import numpy as np
import logging
from . import star#.Star

class Zone():

    """
    Generic class representing a zone.
    """

    def __init__(self,min_radius,max_radius,mass,Z,abunds,*args):
        self.min_r = min_radius
        self.max_r = max_radius
        self.area = np.pi*max_radius**2 - np.pi*min_radius**2
        self.mass = mass
        self.abunds = abunds
        self.Z = Z
        self._stars = []

    def form_star(self,time,sf_mode='sim',mass=mass):
        """
        Adds a new star to the zone.
        """
        
        self._stars.append(Star(time,self.Z,self.abunds,sf_mode=sf_mode,mass))

    def enrich(self,time):
        """
        Loop through all the stars in the zone and enrich all the elements
        we are tracing for one timestep.
        """
        ages = np.array([star.tform - time for star in self._stars])
        iages = np.interp(ages,enrich['times'],np.arange(len(enrich['times'])))
        Zs = np.array([star.Z for star in self._stars])
        iZs = np.interp(Zs,enrich['Zs'],np.arange(len(enrich['Zs'])))
        mstars = np.array([star.mass for star in self._stars])
        for el in self.abunds.keys():
            ej_abunds=scipy.ndimage.map_coordinates(enrich_table[el],zip(iages,iZs),order=1)
            self.abunds[el]=(self.abunds[el]*self.mass + (mstars*ej_abunds).sum()) / self.mass
        rel_masses = scipy.ndimage.map_coordinates(enrich_table['m_ej'],zip(iages,iZs),order=1)
        self.mass += (rel_masses*mstars).sum()
        [star.mass-=star.mass*rel_masses[istar] for istar,star in enumerate(self._stars)]

#        self.abunds,self.Z += snii(tnow, self._stars)
#        self.abunds,self.Z += snia(tnow, self._stars)
#        self.abunds,self.Z += agb(tnow, self._stars)

    def outflow(self):
        """
        Move material out of this zone to target zone
        """

    def inflow(self):
        """
        Add mass from another zone.
        """


def create_zones(n):
    """
    Initializes zone creation, automates creating zones in a disk.
    """
    
    return [Zone(i*dr) for i in range(n)]
