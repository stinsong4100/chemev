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

    def _init_(self,mass,tform,Z,abunds):
        self.mass = mass
        self.tform = tform
        self.Z = Z
        self.abunds = abunds
