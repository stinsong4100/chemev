"""
zones
====

Implements the zone class, which is the fundamental container for gas
in the chemical evolution models.
"""

import numpy as np, scipy.ndimage as spnd
import logging, pdb
from . import enrich, starlifetime as slt, star#.Star

global star_type
star_type = np.dtype({'names':['tform','init_mass','mass','Z'],
                      'formats':['f','f','f','f']})

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

    def form_star(self,time,sf_mode='sim',mass=None):
        """
        Adds a new star to the zone.
        """
        
        new_star = [time,mass,mass,self.Z]
        new_star.extend(self.abunds.values())
        self.stars = np.append(self.stars,np.array([tuple(new_star)],
                                                   dtype=star_type))

    def enrich(self,ages,Zs,s_masses,time_step_length,gas_mass):
        """
        Loop through all the stars in the zone and enrich all the elements
        we are tracing for one timestep.
        """

        max_snii_age = slt.lifetime(8.0,0.02)
        min_snii_age = slt.lifetime(80.,0.02)
        isnii = (ages-time_step_length < max_snii_age) & (ages > min_snii_age)
        iagb = (ages-time_step_length >= max_snii_age)

        self.mass=gas_mass #+=(rel_masses*mstars[isnii]*time_step_length).sum()

#        pdb.set_trace()
        #SNII enrichment
        snii_rel_mass,snii_Zmass = self.enrichment(ages[isnii], Zs[isnii],
                                     s_masses[isnii],'snii',time_step_length)

        #AGB enrichment
        agb_rel_mass,agb_Zmass = self.enrichment(ages[iagb], Zs[iagb],
                                     s_masses[iagb],'agb',time_step_length)

        #SNIa enrichment
        snia_rel_mass,snia_Zmass = self.enrichment(ages[iagb],Zs[iagb],
                                     s_masses[iagb],'snia',time_step_length)

        print "Z Enrichment:  SNII:  %g   AGB:  %g  SNIa:  %g"%(snii_Zmass,agb_Zmass,snia_Zmass)

        TotZmass = self.mass*self.Z + snii_Zmass + agb_Zmass + snia_Zmass
        if (self.mass >0): self.Z = TotZmass / self.mass
        else:  self.Z = 0

        return np.hstack((snii_rel_mass, agb_rel_mass + snia_rel_mass))

    def enrichment(self,ages, Zs,s_masses, fb_type, time_step_length):
        if len(enrich.tables[fb_type]['times'].shape) < 2:
            etimes = enrich.tables[fb_type]['times']
        else:  etimes = enrich.tables[fb_type]['times'][0]

        iages = np.interp(ages,etimes,np.arange(len(etimes)))
        iZs = np.interp(Zs,enrich.tables[fb_type]['Zs'],
                        np.arange(len(enrich.tables[fb_type]['Zs'])))

        for el in self.abunds.keys():
            ej_abund_rates=spnd.map_coordinates(enrich.tables[fb_type]['yield_rates'][el],
                                    np.vstack((iZs,iages)),
                                    order=1,mode='constant',cval=0)
            if self.mass >0: 
                self.abunds[el]=(self.abunds[el]*self.mass + 
                             (s_masses*ej_abund_rates*
                              time_step_length).sum()) / self.mass

        rel_masses = spnd.map_coordinates(enrich.tables[fb_type]['yield_rates']['m_ej'],
                            np.vstack((iZs,iages)),
                            order=1,mode='constant',cval=0)

        Z_masses = spnd.map_coordinates(enrich.tables[fb_type]['yield_rates']['Z'],
                                        np.vstack((iZs,iages)),
                                        order=1,mode='constant',cval=0)
        
        return rel_masses, (s_masses*Z_masses*time_step_length).sum() 

    def outflow(self):
        """
        Move material out of this zone to target zone
        """

    def inflow(self):
        """
        Add mass from another zone.
        """


def create_disk_zones(n,max_disk_r=25,disk_gas_mass=1e10,h_r=4,
                      init_Z=0,init_abunds={'O':0,'Fe':0,'Mg':0}):
    """
    Initializes zone creation, automates creating zones in a disk.
    """

    max_disk_r, h_r = float(max_disk_r), float(h_r)
    disk_zones = []
    dr = max_disk_r / n
    Sigma0 = disk_gas_mass / (np.pi*2*h_r*(h_r - 
                                           np.exp(-max_disk_r/h_r)*
                                           (h_r+max_disk_r)))
    prefac = np.pi*2.0*h_r*Sigma0
    for i in range(n):
        min_r=i*dr
        max_r=(i+1)*dr
        zone_mass = prefac*(np.exp(-min_r/h_r)*(h_r+min_r) - 
                            np.exp(-max_r/h_r)*(h_r+max_r))
    
        disk_zones.append(Zone(min_r,max_r,zone_mass,init_Z,init_abunds))


    return disk_zones
