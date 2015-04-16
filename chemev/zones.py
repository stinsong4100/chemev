"""
zones
====

Implements the zone class, which is the fundamental container for gas
in the chemical evolution models.
"""

import numpy as np, scipy.ndimage as spnd
import logging
from . import data,starlifetime as slt, star#.Star

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
        self.stars = np.array([],dtype=data.star_type)

    def form_star(self,time,sf_mode='sim',mass=None):
        """
        Adds a new star to the zone.
        """
        
        new_star = [time,mass,mass,self.Z].extend(self.abunds.values())
        self.stars = np.append(self.stars,np.array([tuple(new_star)],
                                                   dtype=data.star_type))

    def enrich(self,time,time_step_length,gas_mass):
        """
        Loop through all the stars in the zone and enrich all the elements
        we are tracing for one timestep.
        """
        ages = time - self.stars['tform']
        Zs = self.stars['Z']
        mstars = self.stars['mass']

        max_snii_age = slt.lifetime(8,0.02)
        min_snii_age = slt.lifetime(80,0.02)
        isnii = (ages-time_step_length < max_snii_age) & (ages > min_snii_age)
        iagb = (ages >= max_snii_age)

        #SNII enrichment
        iages = np.interp(ages[isnii],data.snii_dat['times'][0][isort],indices[isort])
        iZs = np.interp(Zs[isnii],data.snii_dat['Zs'][isort],indices[isort])
        Zmass = 0.0
        for el in self.abunds.keys():
            ej_abund_rates=spnd.map_coordinates(data.snii_dat['yield_rates'][el],
                                    np.vstack((iages,iZs)),
                                    order=1,mode='constant',cval=0)
            self.abunds[el]=(self.abunds[el]*self.mass + 
                             (mstars[isnii]*ej_abund_rates*
                              time_step_length).sum()) / self.mass
            if ((el != 'H') or (el !='He')):
                Zmass+=(mstars[isnii]*ej_abund_rates*time_step_length).sum()

        rel_masses = spnd.map_coordinates(data.snii_dat['yield_rates']['m_ej'],
                            np.vstack((iages,iZs)),
                            order=1,mode='constant',cval=0)
        self.stars['mass'][isnii]-=self.stars['mass'][isnii]*rel_masses*time_step_length

        #AGB enrichment
        iages = np.interp(ages[iagb],data.agb_dat['times'][0][isort],indices[isort])
        iZs = np.interp(Zs[iagb],data.agb_dat['Zs'][isort],indices[isort])
        Zmass = 0.0
        for el in self.abunds.keys():
            ej_abund_rates=spnd.map_coordinates(data.agb_dat['yield_rates'][el],
                                    np.vstack((iages,iZs)),
                                    order=1,mode='constant',cval=0)
            self.abunds[el]=(self.abunds[el]*self.mass + 
                             (mstars[iagb]*ej_abund_rates*
                              time_step_length).sum()) / self.mass
            if ((el != 'H') or (el !='He')):
                Zmass+=(mstars[iagb]*ej_abund_rates*time_step_length).sum()

        rel_masses = spnd.map_coordinates(data.agb_dat['yield_rates']['m_ej'],
                            np.vstack((iages,iZs)),
                            order=1,mode='constant',cval=0)
        self.stars['mass'][iagb]-=self.stars['mass'][iagb]*rel_masses*time_step_length


        TotZmass = self.mass*self.Z + Zmass
        self.mass=gas_mass #+= (rel_masses*mstars[isnii]*time_step_length).sum()
        if (self.mass >0): self.Z = TotZmass / self.mass
        else:  self.Z = 0
#        import pdb; pdb.set_trace()

#        self.abunds,self.Z += snii(tnow, self.stars)
#        self.abunds,self.Z += snia(tnow, self.stars)
#        self.abunds,self.Z += agb(tnow, self.stars)

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
