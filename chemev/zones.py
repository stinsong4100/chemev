"""
zones
====

Implements the zone class, which is the fundamental container for gas
in the chemical evolution models.
"""

import numpy as np
import logging

class Zone():

    """
    Generic class representing a zone.
    """

    def __init__(self,radius,mass,Z,*args):
        self.r = radius


    def add_star():
        """
        Adds a new star to the zone.
        """

    def enrich():
        """
        Loop through all the stars in the zone and enrich all the elements
        we are tracing for one timestep.
        """

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
