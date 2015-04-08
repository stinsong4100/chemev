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

    def _init_(self,tform,zonemass,zonearea,Z,abunds):
        self.mass = star_form_efficiency * zonemass#/zonearea
        self.tform = tform
        self.Z = Z
        self.abunds = abunds

