'''
- Physical constants: h, c, k, sigma, G, N_avogadro, Rstart
- Plank function
- Demonstration of how to convert differnt units, using the package *pint*.

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)
'''

import numpy as np
from scipy import constants

# Get basic physical and thermodynamic constants

# Physical constants
h     = constants.h     #Planck's constant
c     = constants.c     #Speed of light
k     = constants.k     #Boltzman thermodynamic constant
sigma = constants.sigma #Stefan-Boltzman constant
G     = constants.G     #Gravitational constant

# Thermodynamic constants

# Following will come out in J/(deg kmol), so that dividing Rstar by molecular
# weight gives gas constant appropriate for mks units
N_avogadro = constants.Avogadro  #Avogadro's number
Rstar = 1000. * k * N_avogadro   #Universal gas constant

#----------------Radiation related functions-------------
def B(nu,T):
    '''Planck function (of frequency)
    Density of blackbody radiation.
    
    Parameters
    ----------
        nu : float
            Frequency [Hz]
        T : float
            Temperature [K]
            
    Return
    ------
        B : float
            Corresponding blackbody radiation density
    
    .. math::
        B(\\nu, T) = \\frac{2 k^3 T^3}{h^2 c^2} \\frac{u^3}{e^u -1}
    
    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>>  
    >>> nu = np.linspace(1e12, 1e14, 1001)
    >>> T = 300
    >>> blackbody_radiation = B(nu, T)
    >>>
    >>> plt.plot(nu, blackbody_radiation)
    >>> plt.xlabel(r'$\\nu [Hz]$')
    >>> plt.ylabel(r'$B(\\nu)$')
    >>> plt.title('Blackbody Radiation')
    >>> plt.show()

    '''

    u = h*nu/(k*T)
    u = np.atleast_1d(u)
    u[u>500] = 500           #To prevent overflow
    B = 2.*h*nu**3/c**2 / (np.exp(u)-1.)
    if len(B)==1:
        B = B[0]
    
    return B

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
    
    
    # Show how to use the Planck-function
    import matplotlib.pyplot as plt
    
    nu = np.linspace(1e12, 1e14, 1001)
    T = 300
    print('B_max at nu= {0:6.3e} Hz'.format(58.8e9*T))
    blackbody_radiation = B(nu, T)
    plt.plot(nu, blackbody_radiation)
    plt.xlabel(r'$\nu [Hz]$')
    plt.ylabel(r'$B(\nu)$')
    plt.title('Blackbody Radiation')
    plt.show()
    
    nu = 1.5e13
    print('The radiation density at {0:5.3e} Hz and {1} K is: {2:8.6e}'.format(nu, T, B(nu, T)))
    
