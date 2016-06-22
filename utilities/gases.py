'''
Physical properties of gases

All units are mks units.
This one replaces the section on gases in the old module "phys.py"

Change Log
----------
06/19/2016 : created, based on "phys.py" [ThH]                

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)

ToDo
----
* Finish putting data into the gas database, perhaps including Van der Waals,
  Shomate and Antoine coefficients, at least in selected cases
* Provide R and Rcp values

'''

import pandas as pd
from io import StringIO

'''
Attribues of gases
------------------
    CriticalPointT :  Critical point temperature [K]
    
    CriticalPointP :  Critical point pressure [Pa]
    
    TriplePointT : Triple point temperature [K]
    
    TriplePointP : Triple point pressure [Pa]
    
    L_vaporization_BoilingPoint : Latent heat of vaporization [J/kg]
        at boiling point.  The so-called "boiling point" is
        the temperature at which the saturation vapor pressure
        equals 1 atmosphere (1.013 bar). For CO2, the "boiling point"
        occurs below the triple point temperature, so the condensed
        phase would not be a liquid. Hence, for CO2 the
        latent heat is given at the arbitrary reference point
        of 253K and 29Pa.
        
    L_vaporization_TriplePoint : Latent heat of vaporization [J/kg]
        at the triple point
        
    L_fusion : Latent heat of fusion [J/kg] at the triple point
    
    L_sublimation : Latent heat of sublimation [J/kg] at triple point
    
    rho_liquid_BoilingPoint : Liquid phase density [kg/m**3]
        at the boiling point
        
    rho_liquid_TriplePoint: Liquid phase density [kg/m**3]
        at the triple point
        
    rho_solid : Solid phase density [kg/m**3] at (or sometimes near)
        the triple point
        
    cp : Gas phase specific heat [J/(kg K)], at 298K and 1 bar
    
    gamma : ratio of specific heat at constant pressure
        to specific heat at constant volume. (Generally
        stated at 298K and 1bar)
        
    MolecularWeight : Molecular weight of the dominant isotope
    
    name : Name of the gas 
    
    formula : Chemical formula (e.g. 'CH4')

    L_vaporization : Default value to use for latent heat of
        vaporization.  Set to triple point value, if available,
        else to boiling point value
        
    rho_liquid : Default value to use for liquid phase density.
        Set to triple point value if available otherwise
        set to boiling point value

    R : Gas constant for the individual gas. Computed from
        other data as Rstar/MolecularWeight, when the update()
        method is called [Currently not given]
        
    Rcp : The adiabatic exponent R/cp. Computed from other
        data when the update() method is called. [Currently not given]
'''

#----------Properties of gases-----------------
#The following are approximate mean values
#for "normal" temperatures and pressures, suitable only
#for rough calculations.

# These data could be also placed into an external file. For simplicity,
# I keep them here, and use "io.StringIO" to treat the long string as a file.

data_string="""
CriticalPointT, CriticalPointP, TriplePointT, TriplePointP, L_vaporization_BoilingPoint, L_vaporization_TriplePoint, L_fusion, L_sublimation, rho_liquid_BoilingPoint, rho_liquid_TriplePoint, rho_solid, cp, gamma, MolecularWeight, name, formula, L_vaporization, rho_liquid
6.4710e2,    2.2100e7,    2.7315e2,    6.1100e2,    2.2550e6,    2.4930e6,    3.3400e5,    2.8400e6,    9.5840e2,    9.9987e2,    9.1700e2,    1.8470e3,    1.3310,    1.8000e1,    Water,          H2O,   2.4930e6,   9.9987e2
1.9044e2,    4.5960e6,    9.0670e1,    1.1700e4,    5.1000e5,    5.3600e5,    5.8680e4,    5.9500e5,    4.5020e2,    None,        5.0930e2,    2.1950e3,    1.3050,    1.6000e1,    Methane,        CH4,   5.3600e5,   4.5020e2
3.0420e2,    7.3825e6,    2.1654e2,    5.1850e5,    None,        3.9700e5,    1.9600e5,    5.9300e5,    1.0320e3,    1.1100e3,    1.5620e3,    8.2000e2,    1.2940,    4.4000e1,    Carbon_Dioxide, CO2,   3.9700e5,   1.1100e3
1.2620e2,    3.4000e6,    6.3140e1,    1.2530e4,    1.9800e5,    2.1800e5,    2.5730e4,    2.4370e5,    8.0860e2,    None,        1.0260e3,    1.0370e3,    1.4030,    2.8000e1,    Nitrogen,       N2,    2.1800e5,   8.0860e2
1.5454e2,    5.0430e6,    5.4300e1,    1.5000e2,    2.1300e5,    2.4200e5,    1.3900e4,    2.5600e5,    1.1410e3,    1.3070e3,    1.3510e3,    9.1600e2,    1.3930,    3.2000e1,    Oxygen,         O2,    2.4200e5,   1.3070e3
3.3200e1,    1.2980e6,    1.3950e1,    7.2000e3,    4.5400e5,    None,        5.8200e4,    None,        7.0970e1,    None,        8.8000e1,    1.4230e4,    1.3840,      2.0000,    Hydrogen,       H2,    4.5400e5,   7.0970e1
5.1000,      2.2800e5,    2.1700,      5.0700e3,    2.0300e4,    None,        None,        None,        1.2496e2,    None,        2.0000e2,    5.1960e3,    1.6640,      4.0000,    Helium,         He,    2.0300e4,   1.2496e2
4.0550e2,    1.1280e7,    1.9540e2,    6.1000e3,    1.3710e6,    1.6580e6,    3.3140e5,    1.9890e6,    6.8200e2,    7.3420e2,    8.2260e2,    2.0600e3,    1.3090,    1.7000e1,    Ammonia,        NH3,   1.6580e6,   7.3420e2
None,        None,        None,        None,        None,        None,        None,        None,        None,        None,        None,        1004.,       1.4003,       28.97,    Earth_Air,      air,   None,       None
"""

gas_props = pd.read_csv(StringIO(data_string), sep='[, ]+', engine='python', na_values='None')

# Allow easy access to the data through the "formula"
gas_props.index = gas_props['formula']

# Create units dictionary
units = ['K', 'Pa', 'K', 'Pa', 'J/kg', 'J/kg', 'J/kg', 'J/kg', 'kg/m**3', 'kg/m**3', 'kg/m**3', 'J/(kg K)', 'None', 'None', 'None', 'None', 'J/kg', 'kg/m**3']
unit_dict = dict(zip(gas_props.columns, units))
    
if __name__ == '__main__':
    # A few examples of how to access the gas-properties
    
    # List of properties stored
    print('The following informations are stored: -------')    
    print(gas_props.columns)    
    
    # Show the properties of Water
    print('The Properties of H2O: -------')
    print(gas_props.loc['H2O'])    
    
    # Show the "cp"-s:
    print('CPs: ---------------')
    print(gas_props.cp)        
        
    # Print all the "Critical*" Parameters    
    print('"Critical" Parameters: ---------------------')
    print(gas_props.filter(regex='Critical*'))    
    
    # Show the units    
    print('Units -----------------')
    for key in unit_dict:
        print('{0}: [{1}]'.format(key, unit_dict[key]))    