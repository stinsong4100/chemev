import scipy.integrate, numpy as np
import totradpres

# Calculation of number and mass of stars that go SN Type Ia in a
# given stellar mass range.  

EPSSNIA = 1e-7
dFracBinSNIa = 0.05
dGamma_2nd = 2.0
dMBmin = 3.0
dMBmax = 16.0
dMSNIImin = 8.0

# Integrand for secondary mass function: probability of a binary
# of mass dMassB having a secondary of MSBinary->dMass2
def dMSIMFSecInt (dMassB, dMass2):
    # A factor of 1/(dMassB*log(10)) is from the definition of IMF as
    # number per log10(M).  Another factor of 1/dMassB is to
    # normalized the integral over mass ratios.
    return (dMass2/dMassB)**dGamma_2nd*totradpres.IMF(dMassB)/dMassB/dMassB/np.log(10.0)

# IMF of secondary in binary system that goes SN Ia
# The distribution of secondary mass ratios is assumed to be a power
# law, mu^dGamma as in Matteucci & Greggio (1986)
def dMSIMFSec(dMass2):
    global dGamma_2nd, dMBin, dFracBinSNIa
    dNorm_2nd = 2.0**(1 + dGamma_2nd)*(1 + dGamma_2nd)  # Normalization
    # of 2ndary 

    Msup = dMass2 + dMSNIImin  # Mass where primary would have gone supernovaII
    # Msup = 16.0;  // Note that the double integral can be tested by
    # uncommenting the above line, setting dFracBinSNIa = 1.0, and
    # comparing NSNIa with the number of stars in the 3.0-16.0 mass interval.
    dMass2_2 = 2.*dMass2  # Minimum mass of binary
    if (dMass2_2 > dMBmin): Minf=dMass2_2
    else:Minf=dMBmin

    dIMFSec = scipy.integrate.romberg(dMSIMFSecInt, Minf, Msup, args=[dMass2],
                                      tol=EPSSNIA)
    dIMFSec *= dFracBinSNIa * dNorm_2nd;
    return dIMFSec;

# calculate number of SN Type Ia a la Raiteri, Villata, Navarro, A&A
# 315, 105, 1996) Returns number of SN Type Ia that occur during
# timestep in which dMassT1 and dMassT2 are masses of stars that end
# their lives at the end and beginning of timestep, respectively
def NSNIa (dMassT1, dMassT2):
    global dMBmax
    # Exclude primaries that go SNII
    if (dMassT1 >= 0.5*dMBmax):
	return 0.0
    if(dMassT2 > 0.5*dMBmax):
	dMassT2 = 0.5*dMBmax
    return scipy.integrate.romberg(dMSIMFSec, dMassT1, dMassT2, tol=EPSSNIA)

