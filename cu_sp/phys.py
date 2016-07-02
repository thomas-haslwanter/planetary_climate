'''
- Physical constants
- Plank function

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)
'''

import numpy as np
from scipy import constants

#-------------Basic physical constants-------------------
h     = constants.h     #Planck's constant
c     = constants.c     #Speed of light
k     = constants.k     #Boltzman thermodynamic constant
sigma = constants.sigma #Stefan-Boltzman constant
G     = constants.G     #Gravitational constant

#-----------Thermodynamic constants------------
# Following will come out in J/(deg kmol), so that dividing Rstar by molecular
# weight gives gas constant appropriate for mks units

N_avogadro = constants.Avogadro  #Avogadro's number
Rstar = 1000.*k*N_avogadro   #Universal gas constant

#----------------Radiation related functions-------------
def B(nu,T):
    '''Planck function (of frequency)
    
    .. math::
        B(\\nu, T) = \\frac{2 k^3 T^3}{h^2 c^2} \\frac{u^3}{e^u -1}
    
    '''

    
    u = min(h*nu/(k*T),500.) #To prevent overflow
    return (2.*h*nu**3/c**2)/(np.exp(u)-1.)    

if __name__=='__main__':
    # Example for working with units -----------------------
    # While the example is for "temperature", most common physical units
    # can be converted.
    import pint
    
    # Select a temperature
    units = pint.UnitRegistry()
    Quantity = units.Quantity
    
    temp_deg = 30                                # Temperature in deg Celsius
    Temperature = Quantity(temp_deg, units.degC)  # Temperature-object
    temp_K = Temperature.to(units.K).magnitude    # Temperature in deg Kelvin
    
    print('{0:5.1f} deg C = {1:5.1f} deg K'.format(temp_deg, temp_K))
