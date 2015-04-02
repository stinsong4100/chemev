"""
starlifetime
====

Implements the starlifetime functions
"""

import numpy as np
import logging

zmin=7e-5
zmax=3e-2
zsol=0.02
a00,a01,a02 =(10.13,0.07547,-0.008084)
a10,a11,a12=(-4.424,-0.7939,-0.1187)
a20,a21,a22=(1.262,0.3385,0.05417)
a0,a1,a2=(0.0,0.0,0.0)

def coefInit(Z):
    Z=np.min([Z,zmax])
    Z=np.max([Z,zmin])
    logZ = np.log10(Z)
    a0 = a00 + a01*logZ + a02*logZ*logZ
    a1 = a10 + a11*logZ + a12*logZ*logZ
    a2 = a20 + a21*logZ + a22*logZ*logZ
    return (a0,a1,a2)

def lifetime(mass,Z):
    a0,a1,a2=coefInit(Z)
    logStarMass = np.log10(mass)
    logLtime = a0 + a1 *logStarMass + a2*logStarMass*logStarMass
    return 10.0**logLtime

def starMass(starLifetime, Z):
    '''Use quadratic equation to find masses corresponding to stellar lifetime'''
#    if(starLifetime <= 0.0): return 1000         #Time can be zero

    a0,a1,a2=coefInit(Z)
    c = a0;
    c -= np.log10(starLifetime);
    b = a1;
    a = a2;
#    if(b*b - 4*a*c < 0.0): return 1000 # time is too small for fitting formula

    logStarMass = (-b - np.sqrt(b*b - 4*a*c))/(2*a);
    logStarMass[np.isnan(logStarMass)]=3
    return 10.**logStarMass;

    
