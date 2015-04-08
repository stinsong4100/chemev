"""
imf
===

Implements various IMF classes
"""

import sys, math, numpy as np, pdb
import getopt, scipy.integrate
from scipy.special import erf

mass=1.0
M_SQRT2 = np.sqrt(2.0)
M_SQRT1_2 = np.sqrt(1.0/2.0)
Bresnan=False
Chabrier=True
mmax = 80

class Chabrier():
    """
    Chabrier IMF
    """
    def __init__(self):
        self.a1=1.923 
        self.sigma=0.69 
        self.mc=0.079 
        self.a2=0.537
        self.b2=-1.3
        self.m2=1.0

    def imf(self,mass):
        if hasattr(mass, "__len__"): mass = np.array(mass)
        if isinstance(mass,np.ndarray):
            dIMF = np.zeros(len(mass))
            im2 = (mass >= self.m2)
            dIMF[im2] = self.a2*mass[im2]**self.b2
            ilm = (mass <= self.m2) & (mass > self.mc)
            dIMF[ilm] = \
                self.a1 * np.exp(- (np.log10(mass[ilm]) - 
                                    np.log10(self.mc))**2.0
                                   /(2.0*self.sigma*self.sigma))
            dIMF[mass < self.mc] = 0.0
            dIMF[mass > mmax] = 0.0
            return dIMF
        else: 
            if (mass > mmax): return 0.0
            elif (mass > self.m2): return self.a2*mass**self.b2;
            elif (mass > self.mc): return self.a1 * np.exp(- (np.log10(mass) - 
                            np.log10(self.mc))**2.0
                            /(2.0*self.sigma*self.sigma))
            else: return 0.0

    def IMFIntlogm(self,logMass):
        mass = 10.0**logMass;
        return mass*self.imf(mass);

    def cum_mass(self,mass):
        if hasattr(mass, "__len__"): mass = np.array(mass)
        if isinstance(mass,np.ndarray):
            dCumM = np.zeros(len(mass))
            xmina=np.zeros(len(mass))
            dCumM[mass > mmax] = 0.0
            im2 = (mass >= self.m2)
            ilm = mass < self.m2
            dCumM[im2] = self.a2/(self.b2+1.0)*(mmax**(self.b2+1.0) - \
                               mass[im2]**(self.b2+1.0))/np.log(10.0)
            dCumM[ilm] = self.a2/(self.b2+1.0)*(mmax**(self.b2+1.0) - \
                                     self.m2**(self.b2+1.0))/np.log(10.0)
            imc = (mass > self.mc) 
            xmina[imc] = np.log10(mass[imc])
            xmina[mass <= self.mc] = np.log10(self.mc)
            dCumM[ilm]+=np.array([scipy.integrate.romberg(self.IMFIntlogm, 
                                                          xmin,np.log10(self.m2)) 
                                  for xmin in xmina[ilm]])
            return dCumM
        else: 
            dCumM = 0.0
            if (mass > mmax): return dCumM
            elif (mass > self.m2):return self.a2/(self.b2+1.0)*(mmax**(self.b2+1.0)-\
                               mass**(self.b2+1.0))/np.log(10.0)
            else:
                if mass < self.mc: xmin=np.log10(self.mc)
                else: xmin = np.log10(mass)

                return self.a2/(self.b2+1.0)*(mmax**(self.b2+1.0)-\
                           self.m2**(self.b2+1.0))/np.log(10.0) + \
                           scipy.integrate.romberg(self.IMFIntlogm, 
                                                   xmin,np.log10(self.m2)) 

    def cum_num(self,mass):
        ymax = (np.log10(self.m2) - np.log10(self.mc))/(M_SQRT2*self.sigma);
        if hasattr(mass, "__len__"): mass = np.array(mass)
        if isinstance(mass,np.ndarray):
            dCumN = np.zeros(len(mass))
            im2 = (mass >= self.m2)
            ilm = mass < self.m2
            dCumN[im2]= self.a2/(self.b2)*(mmax**(self.b2) - 
                                           mass[im2]**(self.b2))/np.log(10.0)
            dCumN[ilm]= self.a2/(self.b2)*(mmax**(self.b2) - 
                                           self.m2**(self.b2))/np.log(10.0)
            ymin = np.zeros(len(mass))
            imc = (mass >= self.mc) 
            ymin[imc] = (np.log10(mass[imc]) - 
                         np.log10(self.mc))/(M_SQRT2*self.sigma);
            ymin[mass < self.mc] = 0.0
            dCumN[ilm] += self.a1*self.sigma*np.sqrt(np.pi)*M_SQRT1_2*(erf(ymax) - 
                                                                    erf(ymin[ilm]))
            return dCumN
        else:
            dCumN = 0.0
            if (mass > mmax): return dCumN
            if (mass > self.m2):return self.a2/(self.b2)*(mmax**(self.b2) - 
                                            mass**(self.b2))/np.log(10.0)
            dCumN = self.a2/(self.b2)*(mmax**(self.b2) - 
                                       self.m2**(self.b2))/np.log(10.0)
#            ymax = (np.log10(self.m2) - np.log10(self.mc))/(M_SQRT2*self.sigma)
            if (mass > self.mc): 
                ymin =(np.log10(mass)-np.log10(self.mc))/(M_SQRT2*self.sigma)
            else: ymin = 0.0
            return dCumN+self.a1*self.sigma*np.sqrt(np.pi)* \
                M_SQRT1_2*(erf(ymax) - erf(ymin))


class PiecewisePowerLaw():
    '''
    Generic piecewise power law IMF.  We allow for 3 pieces.
    By default, it is set to Kroupa 2001 values.  Kroup 2001 only has 2
    pieces, so pieces 2 and 3 are identical.
    parameters from Kroupa 2001, equation 2, and ignoring brown dwarfs,
    Also normalized so that the mass integral is 1. 
    NOTE BENE: Kroupa 2001 has a revised IMF in section 6.2 which is
    different than this; however, below is what is used as the default in
    Starburst99
    (http://www.stsci.edu/science/starburst99/mappings/docs/run.html)
    with the exception that the low mass cutoff is .1 instead of the .08
    below and in the Kroupa paper.
    To convert to the IMF(log10(M)) convention of Miller-Scalo, we
    increase the power law by 1 and multiply the coefficient by
    ln(10.0). See, eg., Chabrier 2003, eq. 2 

    '''

    def __init__(self):
        self.a1=0.22038*2.0*np.log(10.0)
        self.b1=-0.3 
        self.m1=0.08 
        self.a2=0.22038*np.log(10.0)
        self.b2=-1.3
        self.m2=0.5 
        self.a3=0.22038*np.log(10.0)
        self.b3=-1.3
        self.m3=1.0

    def imf(self,mass):
        dIMF[mass > self.m3] = self.a3*mass**self.b3
        dIMF[(mass <= self.m3) & (mass > self.m2)] = self.a2*mass**self.b2
        dIMF[(mass <= self.m2) & (mass > self.m1)] = self.a1*mass**self.b1
        dIMF[mass < self.m1] = 0.0
        return dIMF

    def cum_mass(self,mass):
        #need to add mass < 1 Msun lines
        return self.a3/(self.b3 + 1)*(mmax**(self.b3 + 1) - mass**(self.b3 + 1))

    def cum_num(self,mass):
        # Don't worry about masses below 1 Msun 
        return self.a3/(self.b3)*(mmax**(self.b3) - mass**(self.b3))

    

class Kroupa93(PiecewisePowerLaw):
    """
    Kroupa 1993 IMF
    """
    def __init__(self):
        self.a1=0.3029*1.86606 
        self.b1=-0.3 
        self.m1=0.08 
        self.a2=0.3029
        self.b2=-1.2 
        self.m2=0.5 
        self.a3=0.3029 
        self.b3=-1.7 
        self.m3=1.0
        

def calc_lum(mass,Bresnan=False,Chabrier=False):
    if Chabrier: 
        a1 = 1.923 #0.158 
        sigma = .69
        mc = .079
        a2 = 0.537 #4.43e-2
        b2 = -1.3
        m2 = 1.0
    else:
        a1 = 0.3029*1.86606
        b1 = -0.3
        m1 = .08
        a2 = 0.3029
        b2 = -1.2
        m2 = 0.5
        a3 = 0.3029
        b3 =  -1.7
        m3 = 1.0

    mmax=100

    if Bresnan: 
        La0 = 81.0 
        La1 = 1.78 
        La2 = 0.75
        Lb0 = 2.14
        Lb1 = 3.5
        Lb2 = 4.8
        Lm0 = 20.0
        Lm1 = 2.0
        Lm2 = 1.0
        
        if(mass > Lm0) :
            return a3*La0/(b3 + Lb0)*(mmax**(b3 + Lb0) - mass**(b3 + Lb0))
        else :
            dCumL = a3/(b3 + Lb0)*(mmax**(b3 + Lb0)- Lm0**(b3 + Lb0)) 
        if (mass > Lm1) :
            dCumL += a3*La1/(b3 + Lb1)*(Lm0**(b3 + Lb1) - mass**(b3 + Lb1)) 
            return dCumL
	
        else:
            dCumL += a3*La1/(b3 + Lb1)*(Lm0**(b3 + Lb1) - Lm1**(b3 + Lb1)) 

        if (mass > Lm2) :
            dCumL += a3*La2/(b3 + Lb2)*(Lm1**(b3 + Lb2) - mass**(b3 + Lb2)) 
            return dCumL
        else:
            dCumL += a3*La2/(b3 + Lb2)*(Lm1**(b3 + Lb2) - Lm2**(b3 + Lb2)) 

        if(mass > m3) :
            dCumL += a3*La2/(b3 + Lb2)*(Lm2**(b3 + Lb2) - mass**(b3 + Lb2))
            return dCumL
	
        else:
            dCumL += a3*La2/(b3 + Lb2)*(Lm2**(b3 + Lb2) - m3**(b3 + Lb2))
            return dCumL
        
    else:
    
        La0 = 1.0
        Lpow0 = 4
        La1 = 100
        Lpow1 = 2
        Llim1 = 10
        
        if(mass > mmax):
            return 0
        if Chabrier:
            if (mass > Llim1):
                return a2*La1/(b2 + Lpow1)*(mmax**(b2 + Lpow1)
					  - mass**(b2 + Lpow1))/np.log(10.0)
            else:
                print "a2: "+str(a2)
                dCumL = a2*La1/(b2 + Lpow1)*(mmax**(b2 + Lpow1)
				  - Llim1**(b2 + Lpow1))/np.log(10.0)
                print "M > 10 M_sun dCumL: " + str(dCumL)
            if(mass > m2):
                print "a2: "+str(a2)
                dCumL += a2/(b2 + Lpow0)*(Llim1**(b2 + Lpow0)
				  - mass**(b2 + Lpow0))/np.log(10.0)
            else:
                print "a2: "+str(a2)
                dCumL += a2/(b2 + Lpow0)*(Llim1**(b2 + Lpow0)
				   - m2**(b2 + Lpow0))/np.log(10.0) 
                print "M < 10 M_sun dCumL: " + str(a2/(b2 + Lpow0)*(Llim1**(b2 + Lpow0)
				   - m2**(b2 + Lpow0))/np.log(10.0))
            return dCumL
            # Don't worry about masses below 1 Msun 
        else:
            if(mass > Llim1):
                return a3*La1/(b3 + Lpow1)*(mmax**(b3 + Lpow1)
					- mass**(b3 + Lpow1))
            else:
                dCumL = a3*La1/(b3 + Lpow1)*(mmax**(b3 + Lpow1)
				     - Llim1**(b3 + Lpow1)) 
            if(mass > m3):
                dCumL += a3*La0/(b3 + Lpow0)*(Llim1**(b3 + Lpow0)
					  - mass**(b3 + Lpow0))
            else:
                dCumL += a3*La0/(b3 + Lpow0)*(Llim1**(b3 + Lpow0)
					  - m3**(b3 + Lpow0)) 
                return dCumL

